from django.http import JsonResponse

def addMoleculeJsonResponse(success, data="", error="", smiles=None):
    if success:
        if smiles:
            return JsonResponse({"success": True, "data": data, "smiles": smiles})
        else:
            return JsonResponse({"success": True, "data": data})
    else:
        return JsonResponse({"success": False, "error": error})