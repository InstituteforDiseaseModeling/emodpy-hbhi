import os
from snt.load_paths import load_box_paths

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "..", "..", ".."))

# change to your input folder which contains the required files...
download_dir = os.path.join(PROJECT_DIR, "downloads")
schema_file = os.path.join(download_dir, "schema.json")
eradication_path = os.path.join(download_dir, "Eradication")

# user test data directory
# specify user data for testing
# USER_PATH = r'C:\Projects\emodpy-snt\data'

# load project path
data_path, project_path = load_box_paths(country_name='Example')
input_dir = os.path.join(project_path, "simulation_inputs")

sif_id = '../../../../../dtk_sif.id'
