import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger




input_folder = 'imagejresults/complexes/'
output_folder = 'python_results/colocalisation/'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)


proteins = [folder for folder in os.listdir(f'{input_folder}')]


shsp_trajectories_files=[]
CLIC_trajectories_files=[]
trajectories_files =[[f'{root}/{filename}' for filename in files if 'traj' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
trajectories_files=[item for sublist in trajectories_files for item in sublist]


#this dataframe contains ALL trajectories from both proteins, and the metadata so I can grab specific ones later if i want
molecules = []
for filepath in trajectories_files:


    data = pd.read_csv(filepath)
    data.drop([col for col in data.columns.tolist() if ' ' in col], axis=1, inplace=True)
    data = data.T.reset_index().rename(columns = {'index': 'molecule_number'})


    filepath = filepath.replace('\\', '/')
    file_details = re.split('/|\\\\', filepath)
    timepoint = file_details[-3].split('-')[0]
    colocalisation = file_details[-2]
    protein = file_details[-4]
    
    data['protein'] = protein
    data['timepoint'] = timepoint
    data['colocalisation']=colocalisation

    molecules.append(data)
molecules=pd.concat(molecules)



#now need to find the number of CLIC molecules that are present at each time point
timepoints = molecules['timepoint'].unique().tolist()


shsp=molecules[molecules["protein"] == "hsp27"]
client=molecules[molecules['protein']=='Rhodanese']

#find % colocalisation for shsps over time
zero=[]
twenty=[]
forty=[]
sixty=[]
four_hour=[]
seven_hour=[]
for timepoint in timepoints: 
    timepoint
    protein = 'shsp'
    timepoint_shsp= shsp[shsp["timepoint"] == timepoint]
    total_shsp_timepoint=len(timepoint_shsp)
    coloc_timepoint_shsp=len(timepoint_shsp[timepoint_shsp["colocalisation"] == "coloc"])
    percent_colocal_shsp_timepoint=coloc_timepoint_shsp/total_shsp_timepoint*100
    listo=[timepoint, protein, total_shsp_timepoint,coloc_timepoint_shsp,percent_colocal_shsp_timepoint]
    if 'zero' in timepoint:
        zero.append(listo)
    elif '20min' in timepoint:
        twenty.append(listo)
    elif '40min' in timepoint:
        forty.append(listo)
    elif '60min' in timepoint:
        sixty.append(listo)
    elif '4h' in timepoint:
        four_hour.append(listo)
    elif '7h' in timepoint:
        seven_hour.append(listo)


#find % colocalisation for CLIC over time
for timepoint in timepoints: 
    timepoint
    protein = 'client'
    timepoint_client= client[client["timepoint"] == timepoint]
    total_client_timepoint=len(timepoint_client)
    coloc_timepoint_client=len(timepoint_client[timepoint_client["colocalisation"] == "coloc"])
    percent_colocal_client_timepoint=coloc_timepoint_client/total_client_timepoint*100
    listo=[timepoint, protein, total_client_timepoint,coloc_timepoint_client,percent_colocal_client_timepoint]
    if 'zero' in timepoint:
        zero.append(listo)
    elif '20min' in timepoint:
        twenty.append(listo)
    elif '40min' in timepoint:
        forty.append(listo)
    elif '60min' in timepoint:
        sixty.append(listo)
    elif '4h' in timepoint:
        four_hour.append(listo)
    elif '7h' in timepoint:
        seven_hour.append(listo)




percent_colocalisation_all=[zero,twenty,forty,sixty,four_hour,seven_hour]
column_names=['timepoint', 'protein','total number of proteins','colocalised proteins','percent colocalisation']
percent_colocalisation_all=pd.DataFrame([item for sublist in percent_colocalisation_all for item in sublist])
percent_colocalisation_all.columns=column_names

percent_colocalisation_all.to_csv(f'{output_folder}percent_colocalisation_all.csv')

test=pd.melt(percent_colocalisation_all, id_vars=['timepoint','protein'], value_vars='percent_colocalisation')
ax = sns.barplot(x="timepoint", y="value", hue="protein", data=test, palette='viridis')
ax.set_ylabel('percent_colocalisation')

ax = ax.get_figure()    
ax.savefig(f'{output_folder}percent_colocalisation.png', dpi=400)



