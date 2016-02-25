import moldb.views

"""molinf URL Configuration

The `urlpatterns` list_molecules routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin

urlpatterns = [
    #url(r'^admin/', admin.site.urls),
    url(r'^$', moldb.views.home, name='home'),
    url(r'^add-molecule$', moldb.views.add_molecule, name='add_molecule'),
    url(r'^list-molecules$', moldb.views.list_molecules, name='list_molecules'),

    # API
    url(r'^api/molConverter$', moldb.views.api_molConverter, name='api_molConverter'),
    url(r'^api/addMolecule$', moldb.views.api_addMolecule, name='api_addMolecule')
]