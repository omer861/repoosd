import pandas as pd
from keras.models import load_model
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import pairwise2
import joblib
import os
import nglview as nv
from igfold import IgFoldRunner
from django.conf import settings

# Load data
d = pd.read_csv('Viral.csv')

# Load the model and compile it
model = load_model('viral_antibody.h5')
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mae'])

# Define the extract_features function
def extract_features(sequence):
    analysis = ProteinAnalysis(sequence)
    return analysis.get_amino_acids_percent()

# Process the dataset to extract features for the HC sequences only
features = []
for i in range(len(d)):
    hc_sequence = d['HC'][i]  # Extract HC sequence
    combined_features = extract_features(hc_sequence)
    features.append(combined_features)

feature_df = pd.DataFrame(features).fillna(0)  # Handle any missing values
X = feature_df
# Target/Label
y = d['Pred_affinity']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Scale the features
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Save the scaler for later use
joblib.dump(scaler, 'scaler1.pkl')

def predict_antibody_hc(sequence, X_test, y_test, print_similarity=False):


    # Check similarity for HC sequences
    def calculate_similarity(seq1, seq2):
        alignments = pairwise2.align.globalxx(seq1, seq2)
        max_score = max([alignment[2] for alignment in alignments])
        return max_score / min(len(seq1), len(seq2)) * 100
   
    

    similarities = d['Viral Antigen Sequence'].apply(lambda x: calculate_similarity(sequence, x)) 
    max_similarity = max(similarities)
   
   
    if max_similarity < 85:
        print("Error: The sequence similarity is below 85%.")
        return None, None  # Return a tuple with None values

    if print_similarity:
        print(f"Max similarity score: {max_similarity:.2f}%")

    best_match_index = similarities.idxmax()
    best_match_row = d.iloc[best_match_index]

    hc_sequence = best_match_row['HC']  # Use only the HC sequence for prediction
    new_features = extract_features(hc_sequence)
    new_features_df = pd.DataFrame([new_features]).fillna(0)  # Handle any missing values
    new_features_df = scaler.transform(new_features_df)

   
    predicted_affinity = model.predict(new_features_df)[0][0]

    # Predict on X_test and calculate MAE
    y_pred = model.predict(X_test)
    mae = mean_absolute_error(y_test, y_pred)
    print(f"Mean Absolute Error: {mae:.4f}")

   

    predicted_antibody = {
        'HC': hc_sequence,
        'LC': best_match_row['LC'],
        'CDRH1': best_match_row['CDRH1'],
        'CDRH2': best_match_row['CDRH2'],
        'CDRH3': best_match_row['CDRH3'],
        'CDRL1': best_match_row['CDRL1'],
        'CDRL2': best_match_row['CDRL2'],
        'CDRL3': best_match_row['CDRL3']
    }

    print("Predicted Antibody Sequence:")
    for name, sequence in predicted_antibody.items():
        print(f"{name}: {sequence}")

    print(f"Predicted Affinity for HC: {predicted_affinity:.4f}")

    igfold = IgFoldRunner()

    sequences = {
        "H": f"{predicted_antibody['HC']}{predicted_antibody['CDRH1']}{predicted_antibody['CDRH2']}{predicted_antibody['CDRH3']}",
        "L": f"{predicted_antibody['LC']}{predicted_antibody['CDRL1']}{predicted_antibody['CDRL2']}{predicted_antibody['CDRL3']}"
    }

    try:

        # Define the directory and ensure it exists
        pdb_dir = os.path.join(settings.MEDIA_ROOT, 'pdb_files')
        os.makedirs(pdb_dir, exist_ok=True)

        # Define the full path for the output PDB file
        pdb_file_path = os.path.join(pdb_dir, 'predicted_antibody_structure.pdb')
        igfold.fold(
            pdb_file_path,
            sequences=sequences,
            do_refine=False,
            do_renum=False
        )
        print(f"3D structure prediction completed and saved as '{pdb_file_path}'.")
    except Exception as e:
        print(f"Error in folding: {e}")
        return None, None

    # Check if PDB file was created
    if not os.path.isfile(pdb_file_path):
        print("Error: PDB file not found.")
        return None, None

    return predicted_antibody, predicted_affinity

# Function to visualize the 3D structure using nglview
def visualize_structure_nglview(pdb_file):
    # Create an NGLView widget
    view = nv.show_file(pdb_file)
    view.add_cartoon()
    view.color_by('spectrum')
    return view
