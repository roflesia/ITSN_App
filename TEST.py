

#librerias fijas

import streamlit as st
import pandas as pd
import numpy as np
import PIL as pil

from Bio import PDB
import plotly.express as px


#librerias en prueba


#Page set up-----------------------------------------------------------------------------
st.set_page_config(
    page_title="TEST app",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
    menu_items={
        #pestañas del menu que se pueden modificar, pero solo acepta estas 3
        'Get Help': 'https://www.extremelycoolapp.com/help',
        'Report a bug': "https://www.extremelycoolapp.com/bug",
        'About': "# This is a header. This is an *extremely* cool app!"
    }
)

gene1=("GABPA")


#RUTAS INCOMPLETAS

ruta_af_mut=('ARCHIVOSPDB\AF\MUT')


barra_invertida = '\\'
fin_ruta=('.pdb')

#RUTA COMPLETAs

ruta_completa=ruta_af_mut+barra_invertida+gene1+fin_ruta

Atom_serial_number = []
Atom_name = []
Residue_Name = []
X_orthogonal_coordinates = []
Y_orthogonal_coordinates = []
Z_orthogonal_coordinates = []
B_factor = []


#usando el statement with no es necesario cerrar el archivo
with open(ruta_completa, 'r') as pdb_file:
    #print(pdb_file.read())
    for linea in pdb_file:
        if linea.startswith('ATOM') and linea[13:15] == 'CA':
            Atom_serial_number.append(float(linea[6:11]))
            atomname=str(linea[13:16])
            Atom_name.append(atomname.strip()) #Strip elimina los espacios en blanco en eñ string
            Residue_Name.append(str(linea[17:20]))
            X_orthogonal_coordinates.append(float(linea[30:38]))
            Y_orthogonal_coordinates.append(float(linea[38:46]))
            Z_orthogonal_coordinates.append(float(linea[46:54]))
            B_factor.append(float(linea[60:66]))

pdb_dataframe = pd.DataFrame({'Atom Serial Number': Atom_serial_number,
                              'Atom Name': Atom_name,
                              'Residue Name': Residue_Name,
                              'X orthogonal coordinate': X_orthogonal_coordinates,
                              'Y orthogonal coordinate': Y_orthogonal_coordinates,
                              'Z orthogonal coordinate': Z_orthogonal_coordinates,
                              'B factor': B_factor})

st.dataframe(data=pdb_dataframe)


fig = px.line(pdb_dataframe,
              x='Atom Serial Number',
              y="B factor",
              hover_data=['Residue Name'])
st.plotly_chart(fig, theme="streamlit", use_container_width=True)