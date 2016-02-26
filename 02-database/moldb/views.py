from django.shortcuts import render
from moldb import models
from django.template import RequestContext
from .forms import AddMoleculeForm
from django.http import JsonResponse, HttpResponse
from rdkit import Chem
from .helper_functions import addMoleculeDictSerializer
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

def home(request):
    return render(request,
        'home.html',
        {"molecules_count": models.Molecule.objects.count()})

def add_molecules(request):
    return render(request,
        'add_molecules.html',
        {"form": AddMoleculeForm()})

def list_molecules(request):
    mols_list = models.Molecule.objects.all()
    paginator = Paginator(mols_list, 5)

    page = request.GET.get('page')
    try:
        mols = paginator.page(page)
    except PageNotAnInteger:
        mols = paginator.page(1)
    except EmptyPage:
        mols = paginator.page(paginator.num_pages)

    return render(request,
        'list_molecules.html',
        {"mols": mols, "all_mols_count": len(mols_list)})

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
            return JsonResponse(addMoleculeDictSerializer(False, error="Input format unknown."))

        if mol:
            if format_to == "molfile":
                output = Chem.MolToMolBlock(mol)
            elif format_to == "smiles":
                output = Chem.MolToSmiles(mol)
            else:
                return JsonResponse(addMoleculeDictSerializer(False, error="Output format unknown."))

            if output:
                return JsonResponse({"success": True, "data": output})
        else:
            return JsonResponse(
                addMoleculeDictSerializer(False,
                                          error="Cannot convert from '{}' to '{}': probably invalid data supplied."
                                            .format(format_from, format_to)))

def api_addMolecule(request):
    if request.POST and request.POST["molfile"]:
        return JsonResponse([saveMol(molfile=request.POST["molfile"])], safe=False)

def api_uploadMolecules(request):
    try:
        file = request.FILES["file"]
    except MultiValueDictKeyError:
        return JsonResponse([addMoleculeDictSerializer(False, error="No file supplied.", smiles="-")], safe=False)

    response = []

    for chunk in file.chunks():
        for line in str(chunk, encoding="utf-8").split("\n"):
            response.append(saveMol(smiles=line))

    return JsonResponse(response, safe=False)

def api_downloadMolecules(request):
    folder = "static/temp/"

    if "download_all" in request.GET.keys():
        mols = models.Molecule.objects.all()
        filename = "molecules_all_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))
    else:
        mol_ids = [int(x) for x in request.GET.getlist("mol_ids[]")]
        mols = models.Molecule.objects.filter(id__in=mol_ids)
        filename = "molecules_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))
        #newfile = NamedTemporaryFile(suffix='.sdf')

    path = folder + filename

    writer = AllChem.SDWriter(path)
    for mol in mols:
        writer.write(AllChem.MolFromSmiles(mol.smiles))
    writer.close()

    with open(path, mode="r", encoding="utf-8") as f:
        response = HttpResponse(FileWrapper(f), content_type='application/download')
        response['Content-Disposition'] = 'attachment; filename={}'.format(filename)
        return response

    #wrapper = FileWrapper(newfile)
    #response = HttpResponse(wrapper, mime_type="application/force-download")
    #response['Content-Disposition'] = 'attachment; filename=%s' % os.path.basename(newfile.name)
    #response['Content-Length'] = os.path.getsize(newfile.name)
    #return response

    #p = Popen("rm %s" % filepath, shell=True)"""

    #return JsonResponse({"url": "http://{}/{}".format(request.get_host(), filename)})

# helper functions

def saveMol(smiles=None, molfile=None):
    if smiles or molfile:
        mol = models.Molecule()

        try:
            mol.save(smiles=smiles, molfile=molfile)
        except models.Molecule.MoleculeExistsInDatabase as e:
            return addMoleculeDictSerializer(False, error="Cannot add the molecule: it already exists in database.", smiles=e.smiles)
        except models.Molecule.MoleculeCreationError:
            if smiles:
                return addMoleculeDictSerializer(False, error="Cannot add the molecule: check your structure (valence etc.).", smiles=smiles)
            else:
                return addMoleculeDictSerializer(False, error="Cannot add the molecule: check your structure (valence etc.).")

        return addMoleculeDictSerializer(True, internal_id=mol.internal_id, smiles=mol.smiles)
    else:
        return addMoleculeDictSerializer(False, error="Cannot add empty molecule.")