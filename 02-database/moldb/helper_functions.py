from django.http import JsonResponse

def addMoleculeDictSerializer(success, internal_id=None, error=None, smiles=None):
    if success:
        response = {"success": True, "internal_id": internal_id}
    else:
        response = {"success": False, "error": error}

    if smiles:
        response.update({"smiles": smiles})

    return response
