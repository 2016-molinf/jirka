from django.shortcuts import render
from moldb import models
from django.shortcuts import render_to_response
from django.template import RequestContext
from .forms import AddMoleculeForm
from django.http import JsonResponse
from rdkit import Chem
from .helper_functions import addMoleculeDictSerializer
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
import json
from django.utils.datastructures import MultiValueDictKeyError

def home(request):
    return render_to_response('home.html',
                              {"molecules_count": models.Molecule.objects.count()})

def add_molecule(request):
    return render_to_response('add_molecule.html',
                              {"form": AddMoleculeForm()})

def list_molecules(request):
    mols_list = models.Molecule.objects.all()
    paginator = Paginator(mols_list, 2) # Show 25 contacts per page

    page = request.GET.get('page')
    try:
        mols = paginator.page(page)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        mols = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        mols = paginator.page(paginator.num_pages)

    return render_to_response('list_molecules.html',
                              {"mols": mols})

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