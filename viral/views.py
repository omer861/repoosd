from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse
import os
from . import ml_pred

# Create your views here.
# viral/views.py



def home(request):
    return render(request, 'home.html')

def viral(request):
    return render(request, 'viral.html')

def cancer_view(request):
    return render(request, 'cancer.html')

def res(request):
    antigen_sequence = request.GET.get('ANTIGEN_SEQUENCE', '')
    page_type = request.GET.get('page_type', 'viral')

    context = {}

    try:
        pdb_dir = os.path.join(settings.MEDIA_ROOT, 'pdb_files')
        os.makedirs(pdb_dir, exist_ok=True)

        # Define the full path for the output PDB file
        pdb_file_path = os.path.join(pdb_dir, 'predicted_antibody_structure.pdb')
    except Exception as e:
        print(f"Error in pdb file path: {e}")
        return HttpResponse("Error in pdb file path.", status=404)

    #print('here')

    if antigen_sequence:
        # Call the function to predict the antibody and binding affinity
        predicted_antibody, binding_affinity = ml_pred.predict_antibody_hc(antigen_sequence, ml_pred.X_test, ml_pred.y_test)
        print('predection complete')
        if predicted_antibody is not None and binding_affinity is not None:
            context['predicted_antibody'] = predicted_antibody
            context['binding_affinity'] = binding_affinity

            if os.path.exists(pdb_file_path):
                with open(pdb_file_path, 'r') as pdb_file:
                    pdb_content = pdb_file.read()
                    context['pdb_content'] = pdb_content
            else:
                context['error'] = '3D structure file not found.'

        else:
            context['error'] = 'No antibody sequence found or prediction failed.'
    else:
        context['error'] = 'No antigen sequence provided.'

    return render(request, 'result.html', context)


def display_pdb_3d(request, pdb_filename):
    pdb_file_path = os.path.join(settings.MEDIA_ROOT, 'pdb_files', pdb_filename)

    if os.path.exists(pdb_file_path):
        with open(pdb_file_path, 'r') as pdb_file:
            pdb_content = pdb_file.read()
        return render(request, 'pdb_display_3d.html', {'pdb_content': pdb_content})
    else:
        return HttpResponse("PDB file not found.", status=404)