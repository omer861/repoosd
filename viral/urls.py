# viewer/urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('viral/', views.viral, name='viral'),
    path('cancer/', views.cancer_view, name='cancer'),
    path('res/', views.res, name='res'),
    path('pdb/<str:pdb_filename>/', views.display_pdb_3d, name='display_pdb_3d'),
]