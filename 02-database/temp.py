import os
import django
from moldb.search_paginator import SearchPaginator
from molinf import settings
from django.db.models import Q

os.environ["DJANGO_SETTINGS_MODULE"] = "temp_settings"
django.setup()

from moldb.models import *
from rdkit.Chem import AllChem

"""
ids = [x for x in range(100)]
#ids = [1]
#print(len(ids))
p = SearchPaginator(Molecule, ids, 20)
#print(p.object_ids)
print("paginator_pages_count:", p.pages_count)
print(p.objects_count)
#print(p.get_page(1, return_ids=True))

page = p.get_page(1)
#print(list(page))
print("has_next:", page.has_next)
print("has_previous:", page.has_previous)
print("next_page_number:", page.next_page_number)
print("previous_page_number:", page.previous_page_number)
"""

order = Order.objects.get(id=1)
compounds_modified = []

for order_compound in order.ordercompound_set.all():
    compound = order_compound.compound
    compound.purchased_amount = order_compound.amount
    compounds_modified.append(compound)
    print(order_compound.compound.ordercompound_set.all())

for c in compounds_modified:
    print(c.purchased_amount)