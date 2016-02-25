from django.http import JsonResponse

def addMoleculeJsonResponse(success, data="", error=""):
    if success:
        return JsonResponse({"success": True, "data": data})
    else:
        return JsonResponse({"success": False, "error": error})