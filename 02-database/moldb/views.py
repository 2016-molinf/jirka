# TODO: depot web app, adding/removing amount of compounds, usage statistics etc.
# TODO: IUPAC naming from Chemdoodle web or some Python module

from django.shortcuts import render
from django.template.loader import render_to_string
from moldb import models
from django.template import RequestContext
from .forms import AddMoleculeForm
from django.http import JsonResponse, HttpResponse
from rdkit import Chem
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
import json
from django.utils.datastructures import MultiValueDictKeyError
from django.core.files.temp import NamedTemporaryFile
from wsgiref.util import FileWrapper
from rdkit.Chem import AllChem
import os
import random
from subprocess import Popen
import time
import logging
from django.db.models import Q
from django.db.utils import DataError
from psycopg2 import DataError as DataError_psycopg2
import operator
from functools import reduce
from django_rdkit.models import QMOL, Value
from .settings import *
from .search_paginator import SearchPaginator

logger = logging.getLogger(__name__)

def home(request):
    return render(request,
        'home.html',
        {"molecules_count": models.Molecule.objects.count()})

def add_molecules(request):
    #request.session["uploaded_mols"] = []
    #request.session["last_added_mols_index"] = 0
    #request.session["upload_finished"] = False

    return render(request,
        'add_molecules.html',
        {"form": AddMoleculeForm()})

def list_molecules(request):
    paginator = SearchPaginator(models.Molecule,
                                models.Molecule.objects.all().values_list("id", flat=True),
                                MOLECULES_PER_LIST_PAGE)

    if "page" in request.GET.keys():
        page = paginator.get_page(request.GET["page"])
    else:
        page = paginator.get_page(1)

    return render(request,
        'list_molecules.html',
        {"mols": page, "all_mols_count": paginator.objects_count, "paginator": page, "type": "list"})

"""
def list_molecules(request):
    mols_list = models.Molecule.objects.all()
    paginator = Paginator(mols_list, 20)

    page = request.GET.get('get_page')
    try:
        mols = paginator.page(page)
    except PageNotAnInteger:
        mols = paginator.page(1)
    except EmptyPage:
        mols = paginator.page(paginator.num_pages)

    return render(request,
        'list_molecules.html',
        {"mols": mols, "all_mols_count": len(mols_list), "type": "list"})
"""

def search_molecules(request, *args, **kwargs):
    return render(request,
                  'search_molecules.html',
                  {})

def test(request):
    return render(request,
                  'temp_bootstrap.html',
                  {"mols": models.Molecule.objects.all()})

# API

def api_molConverter(request):
    if request.POST:
        format_from = request.POST["format_from"]
        format_to = request.POST["format_to"]
        data = request.POST["data"]

        if format_from == "molfile":
            mol = Chem.MolFromMolBlock(data)
        elif format_from == "smiles":
            mol = Chem.MolFromSmiles(data)
        else:
            return JsonResponse(addMoleculeDictSerialize(False, error="Input format unknown."))

        if mol:
            if format_to == "molfile":
                output = Chem.MolToMolBlock(mol)
            elif format_to == "smiles":
                output = Chem.MolToSmiles(mol)
            else:
                return JsonResponse(addMoleculeDictSerialize(False, error="Output format unknown."))

            if output:
                return JsonResponse({"success": True, "data": output})
        else:
            return JsonResponse(
                addMoleculeDictSerialize(False,
                                         error="Cannot convert from '{}' to '{}': probably invalid data supplied."
                                         .format(format_from, format_to)))

def api_addMolecule(request):
    if request.POST and request.POST["molfile"]:
        return JsonResponseStatus("success", data=[saveMol(molfile=request.POST["molfile"])])

def api_uploadMolecules(request):
    #request.session["upload_finished"] = False
    #request.session["uploaded_mols"] = []
    #request.session["last_added_mols_index"] = 0

    try:
        filetype = request.POST["filetype"]
    except MultiValueDictKeyError:
        return JsonResponseStatus("error", message="No 'filetype' in POST parameters.")

    try:
        file = request.FILES["file"]
    except MultiValueDictKeyError:
        return JsonResponseStatus("error", message="No file supplied (missing 'file' in FILE parameters).")

    data = []

    try:
        if filetype == "smiles":
            for chunk in file.chunks():
                for line in str(chunk, encoding="utf-8").split("\n"):
                    data.append(saveMol(smiles=line))
                    #request.session["uploaded_mols"].append(saveMol(smiles=line))
        elif filetype == "sdf":
            for mol in AllChem.SDMolSupplier(file.temporary_file_path()):
                data.append(saveMol(rdmol=mol))
                #request.session["uploaded_mols"].append(saveMol(rdmol=mol))
        else:
            return JsonResponseStatus("error", message="Invalid filetype.")
    except Exception as e:
        logger.error(str(e))
        return JsonResponseStatus("error", message="Invalid file supplied.")

    return JsonResponseStatus("success", data=data)

    #request.session["upload_finished"] = True
    #return JsonResponse({"success": True})

def api_uploadMoleculesStatus(request):
    if "last_added_mols_index" in request.session.keys():
        try:
            new_added_mols = request.session["uploaded_mols"][request.session["last_added_mols_index"]:]
            request.session["last_added_mols_index"] = len(new_added_mols)
            response = {"success": True,
                        "mol_number": len(request.session["uploaded_mols"]),
                        "uploaded_mols": new_added_mols}

            if "upload_finished" in request.session.keys():
                if request.session["upload_finished"]:
                    response["upload_finished"] = True
                    request.session["upload_finished"] = False
                    request.session["uploaded_mols"] = []
                    request.session["last_added_mols_index"] = 0

            return JsonResponse(response)
        except IndexError:
            return JsonResponse({"success": False, "error": "IndexError"})

def api_downloadMolecules(request):
    folder = "static/temp/"

    if "download_all" in request.GET.keys():
        mols = models.Molecule.objects.all()
        filename = "molecules_all_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))
    else:
        mol_ids = [int(x) for x in request.GET.getlist("mol_ids[]")]
        mols = models.Molecule.objects.filter(id__in=mol_ids)
        filename = "molecules_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))

    path = folder + filename

    writer = AllChem.SDWriter(path)
    for mol in mols:
        writer.write(AllChem.MolFromSmiles(mol.smiles))
    writer.close()

    with open(path, mode="r", encoding="utf-8") as f:
        response = HttpResponse(FileWrapper(f), content_type='application/download')
        response['Content-Disposition'] = 'attachment; filename={}'.format(filename)
        return response

def api_searchMoleculesByStructure(request):
    if "search_paginators" not in request.session.keys():
        request.session["search_paginators"] = []

    if "page" in request.POST.keys():
        return search_pagination(request, page=request.POST["page"], paginator_number=request.POST["paginator_number"])
    else:
        try:
            mol = Chem.MolFromMolBlock(request.POST["molfile"])

            if mol:
                if request.POST["type"] == "exact":
                    mols_list = models.Molecule.objects.filter(inchi=Chem.MolToInchi(mol))
                elif request.POST["type"] == "substructure":
                    mols_list = models.Molecule.objects.filter(rdmol__hassubstruct=Chem.MolToSmiles(mol))

                paginator = SearchPaginator(models.Molecule,
                                            mols_list.values_list("id", flat=True),
                                            MOLECULES_PER_SEARCH_PAGE)
                request.session["search_paginators"].append({
                    "paginator": paginator,
                    "query": mols_list.query
                })

                return search_pagination(request,
                                         paginator=paginator,
                                         paginator_number=len(request.session["search_paginators"]) - 1)
            else:
                raise models.Molecule.MoleculeCreationError
        except models.Molecule.MoleculeCreationError as e:
            return JsonResponseStatus("error", message=e.message)

"""
def api_searchMoleculesByStructure(request):
    try:
        mol = Chem.MolFromMolBlock(request.POST["molfile"])

        if mol:
            if request.POST["type"] == "exact":
                mols_list = models.Molecule.objects.filter(inchi=Chem.MolToInchi(mol))
            elif request.POST["type"] == "substructure":
                mols_list = models.Molecule.objects.filter(rdmol__hassubstruct=Chem.MolToSmiles(mol))

            if request.POST["form"] == "table":
                mols_count = len(mols_list)

                return JsonResponseStatus("success", data={
                                                        "table": render_to_string(
                                                                    '_molecules_table.html',
                                                                    context={"mols": mols_list, "all_mols_count": mols_count, "type": "search"},
                                                                    request=request),
                                                        "mols_count": mols_count})
        else:
            raise models.Molecule.MoleculeCreationError
    except models.Molecule.MoleculeCreationError as e:
        return JsonResponseStatus("error", message=e.message)
"""

def api_searchMoleculesByFilter(request):
    if "search_paginators" not in request.session.keys():
        request.session["search_paginators"] = []

    if "page" in request.POST.keys():
        return search_pagination(request, page=request.POST["page"], paginator_number=request.POST["paginator_number"])
    else:
        search_fields = ["rdmol", "sum_formula", "inchi", "inchi_key", "smarts", "mw", "mw__lt", "mw__gt"]
        search_type = request.POST["search_type"]

        fields = request.POST.getlist("fields[]")
        values = request.POST.getlist("values[]")
        unknown_fields = []
        search_q = Q()

        if fields and values:
            if "smarts" in fields:
                smarts_index = fields.index("smarts")
                search_q = Q(rdmol__hassubstruct=QMOL(Value(values.pop(smarts_index))))
                del fields[smarts_index]

            for field in fields:
                if field not in search_fields:
                    unknown_fields.append(field)

            if unknown_fields:
                return JsonResponseStatus("error",
                                      message="Unknown search fields: {}. Allowed fields are {}".format(
                                          ["'{}'".format(x) for x in unknown_fields],
                                          ["'{}'".format(y) for y in search_fields]
                                      ))

            predicates = list(zip(fields, values))
            q_list = [Q(x) for x in predicates]

            if search_type == "or":
                if fields and values:
                    search_q = search_q | reduce(operator.or_, q_list)
            elif search_type == "and":
                if fields and values:
                    search_q = search_q & reduce(operator.and_, q_list)

            try:
                mols_list = models.Molecule.objects.filter(search_q)

                paginator = SearchPaginator(models.Molecule,
                                            mols_list.values_list("id", flat=True),
                                            MOLECULES_PER_SEARCH_PAGE)
                request.session["search_paginators"].append({
                    "paginator": paginator,
                    "query": mols_list.query
                })

                return search_pagination(request,
                                         paginator=paginator,
                                         paginator_number=len(request.session["search_paginators"]) - 1)
            except DataError_psycopg2 as e:
                return JsonResponseStatus("error", message=str(e))
            except DataError as e:
                return JsonResponseStatus("error", message=str(e))
        else:
            return JsonResponseStatus("error", message="No search fields supplied.")

"""
def api_searchMoleculesByFilter(request):
    search_fields = ["rdmol", "sum_formula", "inchi", "inchi_key", "smarts", "mw", "mw__lt", "mw__gt"]
    search_type = request.POST["search_type"]

    fields = request.POST.getlist("fields[]")
    values = request.POST.getlist("values[]")
    unknown_fields = []
    search_q = Q()

    if fields and values:
        if "smarts" in fields:
            smarts_index = fields.index("smarts")
            search_q = Q(rdmol__hassubstruct=QMOL(Value(values.pop(smarts_index))))
            del fields[smarts_index]

        for field in fields:
            if field not in search_fields:
                unknown_fields.append(field)

            if unknown_fields:
                return JsonResponseStatus("error",
                                      message="Unknown search fields: {}. Allowed fields are {}".format(
                                          ["'{}'".format(x) for x in unknown_fields],
                                          ["'{}'".format(y) for y in search_fields]
                                      ))

        predicates = list(zip(fields, values))
        q_list = [Q(x) for x in predicates]

        if search_type == "or":
            if fields and values:
                search_q = search_q | reduce(operator.or_, q_list)
        elif search_type == "and":
            if fields and values:
                search_q = search_q & reduce(operator.and_, q_list)

        try:
            mols_list = models.Molecule.__objects_list.filter(search_q)
            mols_count = len(mols_list)
        except DataError_psycopg2 as e:
            return JsonResponseStatus("error", message=str(e))
        except DataError as e:
            return JsonResponseStatus("error", message=str(e))

        if request.POST["form"] == "table":
            return JsonResponseStatus("success", data={
                                                    "table": render_to_string(
                                                                '_molecules_table.html',
                                                                context={"mols": mols_list, "all_mols_count": mols_count, "type": "search"},
                                                                request=request),
                                                    "mols_count": mols_count})
    else:
        return JsonResponseStatus("error", message="No search fields supplied.")
"""

# helper functions

def saveMol(smiles=None, molfile=None, rdmol=None):
    if smiles or molfile or rdmol:
        mol = models.Molecule()

        try:
            mol.save(smiles=smiles, molfile=molfile, rdmol=rdmol)
        except models.Molecule.MoleculeExistsInDatabase as e:
            return addMoleculeDictSerialize(False, error=e.message, smiles=e.smiles)
        except models.Molecule.MoleculeCreationError as e:
            if smiles:
                return addMoleculeDictSerialize(False, error=e.message, smiles=smiles)
            else:
                return addMoleculeDictSerialize(False, error=e.message, smiles="-")

        return addMoleculeDictSerialize(True, internal_id=mol.internal_id, smiles=mol.smiles)
    else:
        return addMoleculeDictSerialize(False, error="Cannot add empty molecule.", smiles="-")

def addMoleculeDictSerialize(success, internal_id=None, error=None, smiles=None):
    if success:
        response = {"success": True, "internal_id": internal_id}
    else:
        response = {"success": False, "error": error}

    if smiles:
        response.update({"smiles": smiles})

    return response

def JsonResponseStatus(status, data=None, message=None):
    if status == "success" and data is not None:
        return JsonResponse({"status": "success", "data": data, "message": None})
    elif status == "error" and message is not None:
        return JsonResponse({"status": "error", "data": None, "message": message})

def search_pagination(request, page=1, paginator=None, paginator_number=0):
    paginator_number = int(paginator_number)

    if not paginator:
        paginator = request.session["search_paginators"][paginator_number]["paginator"]

    try:
        page = paginator.get_page(page)
        return JsonResponseStatus("success", data={
            "table": render_to_string(
                '_molecules_table.html',
                context={"mols": page.objects_list,
                         "type": "search",
                         "paginator": page,
                         "paginator_number": paginator_number},
                request=request),
            "mols_count": paginator.objects_count,
        })
    except SearchPaginator.PageDoesNotExist:
        return JsonResponseStatus("error", message="No molecules found in database.")