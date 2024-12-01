
#-----------------------------------------------------------------------------------------------------------------

#   Librerias utilizadas
import streamlit as st
from streamlit_extras.metric_cards import style_metric_cards

import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import roc_utils as ru
from sklearn import metrics
import scipy.stats as stats #test wilcoxon

import json
import os

import Bio.PDB
from Bio.PDB import PDBParser
from io import StringIO
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select

import urllib.request
from stmol import showmol
import py3Dmol

from http.client import IncompleteRead

#nuevas 12-08
import zipfile
import io   
import requests
from io import BytesIO
from Bio.PDB import PDBParser


# import streamlit as st
# import pandas as pd
# import numpy as np
# from plotly.subplots import make_subplots

# import json
# import os

# import roc_utils as ru
# from sklearn import metrics
# import scipy.stats as stats #test wilcoxon

# import plotly.express as px
# import plotly.graph_objects as go

# from streamlit_extras.metric_cards import style_metric_cards

# import Bio.PDB
# from io import StringIO

# import urllib.request
# from stmol import showmol
# import py3Dmol

# from http.client import IncompleteRead
# import PIL as pil
# import json

# from Bio import PDB
# import plotly.express as px
# import plotly.graph_objects as go
# from plotly import subplots

# from stmol import showmol
# import py3Dmol

#///////////////////////////////////////////////////////////////////////////////////////////////////////

#   FUNCIONES PROPIAS

#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   FUNCION GET_PDB_CA
#Parametros
    #Ruta archivo = string- cualquiera de los genes de la bd
#Return
    #dataframe

# def GET_PDB_CA(ruta):
    
#     Atom_serial_number = []
#     Atom_name = []
#     Residue_Name = []
#     X_orthogonal_coordinates = []
#     Y_orthogonal_coordinates = []
#     Z_orthogonal_coordinates = []
#     B_factor = [] 

#     with open(ruta, 'r') as pdb_file:
#                 for linea in pdb_file:
#                     if linea.startswith('ATOM') and linea[13:15] == 'CA':
#                         Atom_serial_number.append(float(linea[6:11]))
#                         atomname=str(linea[13:16])
#                         Atom_name.append(atomname.strip()) #Strip elimina los espacios en blanco en eñ string
#                         Residue_Name.append(str(linea[17:20]))
#                         X_orthogonal_coordinates.append(float(linea[30:38]))
#                         Y_orthogonal_coordinates.append(float(linea[38:46]))
#                         Z_orthogonal_coordinates.append(float(linea[46:54]))
#                         B_factor.append(float(linea[60:66]))
                                
#                 #DATAFRAME CON LOS VALORES DEL PDB- CA
#                 df_pdb = pd.DataFrame({'Atom Serial Number': Atom_serial_number,
#                                             'Atom Name': Atom_name,
#                                             'Residue Name': Residue_Name,
#                                             'X orthogonal coordinate': X_orthogonal_coordinates,
#                                             'Y orthogonal coordinate': Y_orthogonal_coordinates,
#                                             'Z orthogonal coordinate': Z_orthogonal_coordinates,
#                                             'B factor': B_factor})
#     return df_pdb


def GET_PDB_CA(ruta):
    
    Atom_serial_number = []
    Atom_name = []
    Residue_Name = []
    X_orthogonal_coordinates = []
    Y_orthogonal_coordinates = []
    Z_orthogonal_coordinates = []
    B_factor = [] 


    with urllib.request.urlopen(ruta) as response:
        pdb_file_content = response.read().decode('utf-8').splitlines()
    
    for linea in pdb_file_content:
        if linea.startswith('ATOM') and linea[13:15] == 'CA':
            Atom_serial_number.append(float(linea[6:11]))
            atomname = str(linea[13:16])
            Atom_name.append(atomname.strip())
            Residue_Name.append(str(linea[17:20]))
            X_orthogonal_coordinates.append(float(linea[30:38]))
            Y_orthogonal_coordinates.append(float(linea[38:46]))
            Z_orthogonal_coordinates.append(float(linea[46:54]))
            B_factor.append(float(linea[60:66]))
    
     # Crear un DataFrame con los datos extraídos
    df_pdb = pd.DataFrame({
        'Atom Serial Number': Atom_serial_number,
        'Atom Name': Atom_name,
        'Residue Name': Residue_Name,
        'X orthogonal coordinate': X_orthogonal_coordinates,
        'Y orthogonal coordinate': Y_orthogonal_coordinates,
        'Z orthogonal coordinate': Z_orthogonal_coordinates,
        'B factor': B_factor
    })
    return df_pdb
#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   FUNCION GET_PLDDTS
#Parametros
    #pdb1 = dataframe con el primer pdb
    #ref_1 = string con referencia al pdf anterior
    #pdb2 = dataframe con el segundo pdb
    #ref_2 = string con referencia al pdf anterior
#Return
    #fig_plddt

def GET_PLDDTS(pdb_1, ref_1, pdb_2, ref_2, mutacion_posicion, Sel_pred):
    fig_plddt = go.Figure()

    # Configurar visibilidad inicial según el flag Sel_pred
    visible_pdb_1 = "legendonly" if Sel_pred == "Mutation results" else True
    visible_pdb_2 = "legendonly" if Sel_pred == "Wildtype results" else True
    
    fig_plddt.add_trace(go.Scatter( x=pdb_1.index,
                                    y=pdb_1["B factor"],
                                    name=(ref_1),
                                    line=dict(color='#71b9c7'),
                                    visible=visible_pdb_1))
    fig_plddt.add_trace(go.Scatter(x=pdb_2.index,
                                    y=pdb_2["B factor"],
                                    name=(ref_2),
                                    line=dict(color='#db7093'),
                                    visible=visible_pdb_2))
    
    fig_plddt.add_shape(
        type="line",
        x0=mutacion_posicion,  # Posición de la mutación en el eje x
        y0=0,                  # Empieza desde el mínimo del eje y
        x1=mutacion_posicion,  # Termina en la misma posición en x
        y1=max(pdb_1["B factor"].max(), pdb_2["B factor"].max()),  # Altura máxima en el eje y
        line=dict(color="red", width=2, dash="dash")
    )
    
    fig_plddt.add_trace(go.Scatter(
        x=[mutacion_posicion],
        y=[0],  # Punto ficticio para la leyenda
        mode="lines",
        line=dict(color="red", width=2, dash="dash"),
        name="Mutation's position"
    ))

    fig_plddt.update_layout(
        yaxis_title="PLDDT value",
        xaxis_title="CA atom index" 
    )
    return fig_plddt
#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   FUNCION GET_PAE_DF
#Parametros
    #ruta al archivo .json
#Return
    #data frame con las distancias 
# def GET_PAE_DF(ruta):
#     with open(ruta, 'r') as json_file:
#             file_mut=json_file
#             data_mut=json.load(file_mut)
#     df_pae=(pd.DataFrame(data_mut["distance"]))
#     return df_pae


# def GET_PAE_DF(url):
#     # Abrir la URL y leer el contenido
#     with urllib.request.urlopen(url) as response:
#         # Decodificar el contenido y cargarlo como JSON
#         data_mut = json.loads(response.read().decode('utf-8'))
    
#     # Convertir la sección "distance" en un DataFrame
#     df_pae = pd.DataFrame(data_mut["distance"])
#     return df_pae


def GET_PAE_DF(url):
    try:
        # Intentar abrir la URL y leer el contenido
        with urllib.request.urlopen(url) as response:
            # Decodificar el contenido y cargarlo como JSON
            data_mut = json.loads(response.read().decode('utf-8'))
        
        # Convertir la sección "distance" en un DataFrame
        df_pae = pd.DataFrame(data_mut["distance"])
        return df_pae

    except:
        df_pae = pd.DataFrame()
        return df_pae
#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   FUNCION GET_PAE_GRAPH
#Parametros
#   df_pae: dataframe con las distancias del pae
#Return
#   objeto fig_pae 

def GET_PAE_GRAPH(df_pae, mutacion_posicion):
    fig_pae = go.Figure(data=go.Heatmap(
                    z=df_pae, 
                    colorscale="Rainbow"))
    # Obtener la cantidad de filas y columnas del dataframe
    # num_rows, num_cols = df_pae.shape
    tick_values = [mutacion_posicion]  # Puedes ajustar según las posiciones que te interesen
    tick_text = [str(mutacion_posicion)]
    # Ajustar el layout para hacer el heatmap cuadrado
    fig_pae.update_layout(
        width=500,  # Establecer el ancho deseado
        height=500,  # Establecer la altura deseada
        yaxis=dict(scaleanchor="x", scaleratio=1, 
        tickmode='array',
        tickvals=tick_values,
        ticktext=tick_text),
        xaxis=dict(
        tickmode='array',
        tickvals=tick_values,
        ticktext=tick_text
    )
        
    )
    return fig_pae

def GET_PAE_GRAPH_json(distance_matrix):
    fig_pae = go.Figure(data=go.Heatmap(
        z=distance_matrix,
        colorscale="Rainbow"
    ))
    
    # Ajustar el layout para hacer el heatmap cuadrado
    fig_pae.update_layout(
        width=500,  # Establecer el ancho deseado
        height=500,  # Establecer la altura deseada
        yaxis=dict(scaleanchor="x", scaleratio=1),
        xaxis=dict(constrain='domain')
    )
    return fig_pae


#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   FUNCION GET_PAE_DIF
#Parametros 
#   df_1 
#   df_2
#return
#   objeto fig_pae_dif

def GET_PAE_DIF(df_1, df_2):
    dif_pae = df_1 - df_2
    fig_pae_dif = go.Figure(data=go.Heatmap(z=dif_pae, 
                                            colorscale="RdBu")) 
    fig_pae_dif.update_layout(
        width=500,  # Establecer el ancho deseado
        height=500,  # Establecer la altura deseada
        yaxis=dict(scaleanchor="x", scaleratio=1),
        xaxis=dict(constrain='domain')
        
    )
    return fig_pae_dif
#   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#///////////////////////////////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////////////////////////////
# def superimpose_proteins(p1, p2, silent=False):
#     # Protein superposition
#     # Thanks to Anders Steen Christensen for the code: https://gist.github.com/andersx/6354971
#     pdb_parser = Bio.PDB.PDBParser()
#     ref_structure = pdb_parser.get_structure("reference", p1) 
#     sample_structure = pdb_parser.get_structure("sample", p2)

#     # In case personalized results are loaded and there are various models, only the first is chosen
#     # change the value manually if needed!
#     ref_model = ref_structure[0]
#     sample_model = sample_structure[0]

#     ref_atoms = []
#     sample_atoms = []

#     for ref_chain in ref_model:
#         for ref_res in ref_chain:
#             ref_atoms.append(ref_res['CA'])

#     for sample_chain in sample_model:
#         for sample_res in sample_chain:
#             sample_atoms.append(sample_res['CA'])

#     # numpy works better so transforming
#     ref_atoms = np.array(ref_atoms)
#     sample_atoms = np.array(sample_atoms)

#     super_imposer = Bio.PDB.Superimposer()
#     super_imposer.set_atoms(ref_atoms, sample_atoms)
#     super_imposer.apply(sample_structure.get_atoms())  # modifies th

#     # Save the superimposed protein to a string
#     io = Bio.PDB.PDBIO()
#     io.set_structure(sample_structure)
#     output_str = StringIO()
#     io.save(output_str)
#     pdb_content = output_str.getvalue()

#     return pdb_content


def superimpose_proteins(wt_pdb, mut_pdb):
    pdb_parser = Bio.PDB.PDBParser()
    ref_structure = pdb_parser.get_structure("reference", StringIO(wt_pdb))
    sample_structure = pdb_parser.get_structure("sample", StringIO(mut_pdb))

    ref_model = ref_structure[0]
    sample_model = sample_structure[0]

    ref_atoms = [atom for atom in ref_model.get_atoms() if atom.get_id() == 'CA']
    sample_atoms = [atom for atom in sample_model.get_atoms() if atom.get_id() == 'CA']

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_structure.get_atoms())

    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    output_str = StringIO()
    io.save(output_str)
    pdb_content = output_str.getvalue()

    return pdb_content


def superimpose_pdb_with_threshold(wt_url, pred_url):
    """
    Superpone dos estructuras PDB basadas en átomos CA con un pLDDT mayor a 'threshold' 
    y muestra la estructura alineada en Streamlit.

    :param wt_url: URL del archivo PDB de la estructura experimental.
    :param pred_url: URL del archivo PDB de la estructura predicha.
    :param threshold: Umbral de pLDDT para seleccionar residuos.
    :return: String con el contenido PDB de la estructura superimpuesta.
    """
    threshold = 70
    # Cargar las estructuras PDB desde las URLs
    parser = PDBParser(QUIET=True)
    with urllib.request.urlopen(wt_url) as response:
        wt_pdb_data = response.read().decode('utf-8')
        structure_exp = parser.get_structure('exp', StringIO(wt_pdb_data))

    with urllib.request.urlopen(pred_url) as response:
        pred_pdb_data = response.read().decode('utf-8')
        structure_pred = parser.get_structure('pred', StringIO(pred_pdb_data))

    # Obtener átomos CA con el umbral de pLDDT
    exp_atoms, pred_atoms = [], []
    for pred_chain in structure_pred.get_chains():
        for pred_res in pred_chain.get_residues():
            if pred_res.has_id('CA') and pred_res['CA'].get_bfactor() >= threshold:
                pred_res_id = pred_res.get_id()[1]
                for exp_chain in structure_exp.get_chains():
                    for exp_res in exp_chain.get_residues():
                        exp_res_id = exp_res.get_id()[1]
                        if exp_res_id == pred_res_id and exp_res.has_id('CA'):
                            exp_atoms.append(exp_res['CA'])
                            pred_atoms.append(pred_res['CA'])
                            break

    # Realizar la superposición
    super_imposer = Superimposer()
    super_imposer.set_atoms(exp_atoms, pred_atoms)
    super_imposer.apply(structure_pred.get_atoms())

    # Guardar la estructura en un StringIO y devolver su contenido
    io = PDBIO()
    io.set_structure(structure_pred)
    output_str = StringIO()
    io.save(output_str)
    pdb_content = output_str.getvalue()

    return pdb_content
#///////////////////////////////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////////////////////////////
def aadistance(reference_pdb, sample_pdb, mutacion_posicion):
    """
    Vamos a calcular las distancias euclidianas entre los átomos CA de dos estructuras PDB.
    
    Parámetros:
    - reference_pdb (str): Ruta al archivo PDB de referencia.
    - sample_pdb (str): Ruta al archivo PDB a calcular distancias, se sugiere que vaya la superimposed aca! 
    
    Retorna:
    - np.ndarray: Matriz de distancias euclidianas entre los átomos CA de las dos estructuras.
    """
    pdb_parser = PDBParser()

    # Cargar las estructuras PDB
    # ref_structure = pdb_parser.get_structure("reference", reference_pdb)
    # sample_structure = pdb_parser.get_structure("sample", sample_pdb)

    ref_structure = pdb_parser.get_structure("reference", StringIO(reference_pdb))
    sample_structure = pdb_parser.get_structure("sample", StringIO(sample_pdb))

    # Obtener los átomos CA de las estructuras
    ref_atoms = [atom for atom in ref_structure.get_atoms() if atom.get_id() == 'CA']
    sample_atoms = [atom for atom in sample_structure.get_atoms() if atom.get_id() == 'CA']
    
    if len(ref_atoms) != len(sample_atoms):
        raise ValueError("El número de átomos CA en las estructuras no coincide.")
    
    # Obtener las coordenadas de los átomos CA
    ref_coords = np.array([atom.get_coord() for atom in ref_atoms])
    sample_coords = np.array([atom.get_coord() for atom in sample_atoms])
    
    # Calcular las distancias euclidianas
    distancia_euclidiana = np.sqrt(np.sum((ref_coords - sample_coords) ** 2, axis=1))

    fig = go.Figure(data=[go.Bar(
        x=list(range(len(distancia_euclidiana))), 
        y=distancia_euclidiana, 
        marker=dict(color='#71B9C7'),
        name="Euclidean distances"  # Leyenda para las distancias euclidianas
    )])
        # Añadir una línea vertical para la posición de la mutación
    fig.add_shape(
        type="line",
        x0=mutacion_posicion,  # Posición de la mutación en el eje x
        y0=0,                  # Empieza desde el mínimo del eje y
        x1=mutacion_posicion,   # Termina en la misma posición en x
        y1=max(distancia_euclidiana),  # Altura máxima del gráfico
        line=dict(color="red", width=2, dash="dash")  # Línea roja discontinua
    )

    
    fig.add_trace(go.Scatter(
        x=[mutacion_posicion],
        y=[0],  # Punto ficticio para la leyenda
        mode="lines",
        line=dict(color="red", width=2, dash="dash"),
        name="Mutation position"
    ))

   
    fig.update_layout(title="Euclidean distances between CA atoms",
                      xaxis_title="CA atom index",
                      yaxis_title="Euclidean distance",
                      showlegend=True)
    
    
    # Calcular TM-score y RMSD global
    d0_ltarget = 1.24 * np.cbrt(len(ref_atoms) - 15) - 1.8
    tm_score = np.sum(1 / (1 + (distancia_euclidiana / d0_ltarget) ** 2)) / len(ref_atoms)
    rmsd_global = distancia_euclidiana.mean()

    return fig, tm_score, rmsd_global
#///////////////////////////////////////////////////////////////////////////////////////////////////////

def encontrar_posicion_mutacion(secuencia_completa, neoantigeno, posicion_mutada):
    """
    Encuentra la posición de la mutación en la secuencia completa dado un neoantígeno y la posición de la mutación dentro de él.

    Args:
        secuencia_completa (str): Secuencia completa de la proteína.
        neoantigeno (str): Secuencia del neoantígeno (9-10 aminoácidos).
        posicion_mutada (int): Posición del aminoácido mutado dentro del neoantígeno (1-indexed).

    Returns:
        int: Posición de la mutación en la secuencia completa (1-indexed).
        None: Si el neoantígeno no se encuentra en la secuencia completa.
    """
    # Encontrar el índice del neoantígeno en la secuencia completa (0-indexed)
    indice_neoantigeno = secuencia_completa.find(neoantigeno)
    
    if indice_neoantigeno == -1:
        # Si no se encuentra el neoantígeno
        return None

    # La posición de la mutación en la secuencia completa (1-indexed)
    posicion_mutacion = indice_neoantigeno + posicion_mutada

    return posicion_mutacion







#-----------------------------------------------------------------------------------------------------------------
#   Lectura de la base de datos
df_ITSNdb_complete =(pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/ITSNdb.csv", index_col=False))

HLA_df = df_ITSNdb_complete['HLA']
HLAGENE = []
HLAALLELO = []
HLAPROT = []
for n in HLA_df:
    #st.write(n)
    Lc = [*n]
    #st.write(Lc)
    #st.write(Lc[4])
    HLAGENE.append(Lc[4])
    HLAALLELO.append(Lc[5]+Lc[6])
    HLAPROT.append(Lc[8]+Lc[9])

df_ITSNdb_complete["HLA_GENE"] = HLAGENE
df_ITSNdb_complete["HLA_ALLELE_GROUP"] = HLAALLELO
df_ITSNdb_complete["HLA_PROTEIN"] = HLAPROT

#-----------------------------------------------------------------------------------------------------------------
#   Lectura del archivo con info adicional de genes - ITSNdb_GENES.csv

df_ITSNdb_GENES = pd.read_csv("ITSNdb_GENES.csv")

#-----------------------------------------------------------------------------------------------------------------
#   Set Up de la applicación
st.set_page_config(
    page_title="ITSNdb app",
    page_icon="Logo_cel.png",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        # #pestañas del menu que se pueden modificar, pero solo acepta estas 3
        # 'Get Help': "COMPLETAR MAIL",
        # 'Report a bug': "COMPLETAR - MAIL",
        # 'About': "COMPLETAR"
    }
)

#Logo
st.image("logo.png",
         width = 300)


#-----------------------------------------------------------------------------------------------------------------
# SideBar con logos

with st.sidebar:
    st.image("logo.png",
         width = 150)
    st.image("logo_fcefyn.png",
         width = 150)
    st.image("UNC_logo.png",
         width = 150)
    st.image("Logo_conicet.png",
         width = 150)
    st.image("LOGO DATALAB OK.png",
         width = 150)
    st.image("Copia de FPM_Logotipos_RGB-01.png",
         width = 150) 
    
    st.subheader("Usefull links")

    # Botón para el artículo de investigación
    st.page_link(page = "https://doi.org/10.3389/fimmu.2023.1094236", label = "Research article", icon="📑")


    # Botón para el repositorio de GitHub
    st.page_link(page = "https://github.com/roflesia/ITSN_App", label = "GitHub repository", icon="🔗")

#   TABS
tab1, tab2, tab3= st.tabs(["HOME", "DATABASE", "SOFTWARE VALIDATION"])

#-----------------------------------------------------------------------------------------------------------------
#   HOME

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   Previos 
columna = df_ITSNdb_complete.columns 
    #Total entries
Total_entries = df_ITSNdb_complete.shape[0]
    #Tumor Types
Tumor_types = df_ITSNdb_complete.value_counts(columna[0]).shape[0]
    #Positives
Positive_entries = df_ITSNdb_complete.value_counts(columna[6]).iat[0]
    #Negatives
Negative_entries = df_ITSNdb_complete.value_counts(columna[6]).iat[1]
    #Anchor
Anchor_entries = df_ITSNdb_complete.value_counts(columna[8]).iat[1]
    #Non- Anchor
Non_Anchor_entries = df_ITSNdb_complete.value_counts(columna[8]).iat[0]


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   Diseño

with tab1:
    col1, col2 = st.columns(2)
    
    #Columna 1: Texto
    with col1:
        st.header("What is ITSNdb?")

        # Botón para el artículo de investigación
        st.page_link(page = "https://doi.org/10.3389/fimmu.2023.1094236", label = "Access the original research article", icon="📑")
        # Botón para el repositorio de GitHub
        st.page_link(page = "https://github.com/roflesia/ITSN_App", label = "Visit the GitHub repository for this project", icon="🔗")
        
        multi = '''The ITSNdb is a *new* neoantigen database with know immunogenic and non immunogenic tumor specific antigenic peptides derived from genomic rearrangements, such as single nucleotide variants (SNVs), that satisfy the following criteria: '''
        multi2 = '''
          1. The wild type counterpart has been identified in the source protein
          2. The MHC-I presentation has been experimentally validated
          3. The positive or negative immunogenicity has been experimentally validated by, for instance, ELISPOT®'''
        
        multi3 = '''In this sence, all peptides in the database have experimental confirmation of their positive/negative immunogenicity (classified as “Positive” and “Negative” neoantigens respectively) as well as their cell surface presentation.The neoantigens were collected and curated from published articles searched on PubMedTM using “neoantigen'' or “neoepitopes” as keywords. The ITSNdb only includes neoantigens whose inclusion criteria were explicitly described in its reference bibliography. '''
        
        st.markdown(multi, unsafe_allow_html= True)
        st.markdown(multi2, unsafe_allow_html= True)
        st.markdown(multi3, unsafe_allow_html= True)


    #Columna 2: metricas    
    with col2: 
        st.subheader("About the data")
        col3, col4 = st.columns(2)
        with col3:
            st.metric(label = ("Total number of neoantigens:"),
                      value = Total_entries)
            st.metric(label = ("Number of positive neoantigens:"),
                      value = Positive_entries)
            st.metric(label = ("Number of anchor mutations:"),
                      value = Anchor_entries)
        with col4:
            st.metric(label = ("Tumor types:"),
                      value = Tumor_types)
            st.metric(label = ("Number of negative neoantigens:"),
                      value = Negative_entries)
            st.metric(label = ("Number of non-anchor mutations:"),
                      value = Non_Anchor_entries)
        style_metric_cards()

        st.write(" ")


    #Gráficos:
    # grafico de torta con tipos tumorales
    dftumor = df_ITSNdb_complete.groupby(['Tumor']).size().reset_index('Tumor')
    
    fig_tumor = px.pie(dftumor, 
                       names = 'Tumor', 
                       values = 0)
    
    # grafico representacion hla positive/negative
    HLA_NeoType = df_ITSNdb_complete.groupby(['HLA', 'NeoType']).size()
    HLA_NeoType = HLA_NeoType.reset_index('NeoType')
    HLA_NeoType = HLA_NeoType.reset_index('HLA')
    HLA_NeoType = HLA_NeoType.rename(columns={0:' '})

    fig_HLA = px.bar(HLA_NeoType,
                     x="HLA",
                     y=" ",
                     color="NeoType",
                     barmode='relative')

    # grafico representacion hla positive/negative
    PositionType_NeoType = df_ITSNdb_complete.groupby(['PositionType', 'NeoType']).size()
    PositionType_NeoType = PositionType_NeoType.reset_index('PositionType')
    PositionType_NeoType = PositionType_NeoType.reset_index('NeoType')
    PositionType_NeoType = PositionType_NeoType.rename(columns={0:' '}) 
   
    fig_PositionType_NeoType = px.bar(PositionType_NeoType, 
                                      x = "PositionType", 
                                      y = " ", 
                                      color="NeoType",
                                      barmode='relative')
   
    # grafico representacion hla positive/negative
    mutPosition_NeoType = df_ITSNdb_complete.groupby(['mutPosition', 'NeoType']).size()
    mutPosition_NeoType = mutPosition_NeoType.reset_index('mutPosition')
    mutPosition_NeoType = mutPosition_NeoType.reset_index('NeoType')
    mutPosition_NeoType = mutPosition_NeoType.rename(columns={0:' '}) 
   
    fig_mutPosition_NeoType = px.bar(mutPosition_NeoType, 
                                      x = "mutPosition", 
                                      y = " ", 
                                      color ="NeoType",
                                      barmode ='relative')
    
    # grafico representacion hla positive/negative
    Length_NeoType = df_ITSNdb_complete.groupby(['Length', 'NeoType']).size()
    Length_NeoType = Length_NeoType.reset_index('Length')
    Length_NeoType = Length_NeoType.reset_index('NeoType')
    Length_NeoType = Length_NeoType.rename(columns={0:' '}) 
   
    fig_Length_NeoType = px.bar(Length_NeoType, 
                                      x = "Length", 
                                      y = " ", 
                                      color ="NeoType",
                                      barmode ='relative')

    # matriz long 9
        #positivos
        #negativos 
    # matriz long 10
        #positivos
        #negativos

    # Codigo matrices de frecuencia
    #filtrar datos segun longitud e inmunogenicidad

    df_itsn_9 = df_ITSNdb_complete[df_ITSNdb_complete["Length"] == 9]
    df_itsn_9_positives = df_itsn_9[df_itsn_9["NeoType"] == "Positive"]
    df_itsn_9_negatives = df_itsn_9[df_itsn_9["NeoType"] == "Negative"]

    df_itsn_10 = df_ITSNdb_complete[df_ITSNdb_complete["Length"] == 10]
    df_itsn_10_positives = df_itsn_10[df_itsn_10["NeoType"] == "Positive"]
    df_itsn_10_negatives = df_itsn_10[df_itsn_10["NeoType"] == "Negative"]

    #crear listas con los neoantigenos de cada longitud

    list_itsn_9_positives = df_itsn_9_positives["Neoantigen"].to_list()
    list_itsn_9_negatives = df_itsn_9_negatives["Neoantigen"].to_list()
    list_itsn_10_positives = df_itsn_10_positives["Neoantigen"].to_list()
    list_itsn_10_negatives = df_itsn_10_negatives["Neoantigen"].to_list()


    #Lista aminoacidos
    AAm = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    #Crear diccionarios

    dict_9_positives = {}
    dict_9_negatives = {}
    dict_10_positives = {}
    dict_10_negatives = {}

    #Rellenar s/ frecuencias
    #Diccionario de neoantigenos de longitud 9 e inmunogenicos
    for key_pos in range(1,10):
        dict_9_positives[(key_pos)] = {}
        for key_aa in AAm: 
            dict_9_positives[(key_pos)][key_aa] = 0
    
    for i in list_itsn_9_positives:
        for n in range(len(i)): 
            dict_9_positives[(n+1)][i[n]] += 1

    #Diccionario de neoantigenos de longitud 9 y no inmunogenicos
    for key_pos in range(1,10):
        dict_9_negatives[(key_pos)] = {}
        for key_aa in AAm: 
            dict_9_negatives[(key_pos)][key_aa] = 0
    
    for i in list_itsn_9_negatives:
        for n in range(len(i)): 
            dict_9_negatives[(n+1)][i[n]] += 1

    #Diccionario de neoantigenos de longitud 10 e inmunogenicos
    for key_pos in range(1,11):
        dict_10_positives[(key_pos)] = {}
        for key_aa in AAm: 
            dict_10_positives[(key_pos)][key_aa] = 0
    
    for i in list_itsn_10_positives:
        for n in range(len(i)): 
            dict_10_positives[(n+1)][i[n]] += 1

    #Diccionario de neoantigenos de longitud 10 y no inmunogenicos
    for key_pos in range(1,11):
        dict_10_negatives[key_pos] = {}
        for key_aa in AAm: 
            dict_10_negatives[key_pos][key_aa] = 0
    
    for i in list_itsn_10_negatives:
        for n in range(len(i)): 
            dict_10_negatives[(n+1)][i[n]] += 1

    # Dataframes y arrays de frecuencia
    df_9_pos_frec = pd.DataFrame.from_dict(dict_9_positives)
    array_9_pos_frec = df_9_pos_frec.to_numpy()

    df_9_neg_frec = pd.DataFrame.from_dict(dict_9_negatives)
    array_9_neg_frec = df_9_neg_frec.to_numpy()

    df_10_pos_frec = pd.DataFrame.from_dict(dict_10_positives)
    array_10_pos_frec = df_10_pos_frec.to_numpy()

    df_10_neg_frec = pd.DataFrame.from_dict(dict_10_negatives)
    array_10_neg_frec = df_10_neg_frec.to_numpy()

    #arrays de diferencias
    df_dif_9 = df_9_pos_frec - df_9_neg_frec
    array_dif_9 = np.subtract(array_9_pos_frec, array_9_neg_frec)
    array_dif_10 = np.subtract(array_10_pos_frec, array_10_neg_frec)

    #Concatenacion de arrays

    array_9_concat = np.concatenate( (array_9_pos_frec, array_9_neg_frec, array_dif_9), axis = 1)

    vmin = np.amin(array_9_concat)
    vmax = np.amax(array_9_concat)

    #creacion plot 9

    fig_mat_frec_9 = make_subplots(rows = 1, cols = 3, shared_yaxes = True, subplot_titles = ("Inmunogenic", "Non inmunogenic", "Difference"))
    fig_mat_frec_9.add_trace(
    go.Heatmap(z = array_9_pos_frec, 
                x = ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
                y = AAm, 
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=1)
    fig_mat_frec_9.add_trace(
    go.Heatmap(z = array_9_neg_frec, 
                x = ["1", "2", "3", "4", "5", "6", "7", "8", "9"],
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=2)
    fig_mat_frec_9.add_trace(
    go.Heatmap(z = array_dif_9, 
                x = ["1", "2", "3", "4", "5", "6", "7", "8", "9"], 
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=3)

    fig_mat_frec_9.update_layout(coloraxis={'colorscale':'rdbu', 'cmid': 0})

    # plot 10

    fig_mat_frec_10 = make_subplots(rows = 1, cols = 3, shared_yaxes = True, subplot_titles = ("Inmunogenic", "Non inmunogenic", "Difference"))
    fig_mat_frec_10.add_trace(
    go.Heatmap(z = array_10_pos_frec, 
                x = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'],
                y = AAm, 
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=1)
    fig_mat_frec_10.add_trace(
    go.Heatmap(z = array_10_neg_frec, 
                x = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'],
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=2)
    fig_mat_frec_10.add_trace(
    go.Heatmap(z = array_dif_10, 
                x = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'],
                #colorscale = "rdbu",
                coloraxis='coloraxis'),
    row=1, col=3)

    fig_mat_frec_10.update_layout(coloraxis={'colorscale':'rdbu', 'cmid': 0})

    
    if "contador" not in st.session_state:
        st.session_state["contador"] = 0
    if "totalgraficos" not in st.session_state:
        st.session_state["totalgraficos"] = 7


    ColA, ColB, ColC= st.columns([1,4,1], gap="large")

    with ColA:  
        st.header("Statistical plots")
        ColAA, ColAB = st.columns(2)
        with ColAA:
            previous = st.button(":rewind:", use_container_width=True)
        with ColAB:
            next = st.button(":fast_forward:", use_container_width=True)
    
         
    if next:
        st.session_state["contador"] += 1
    if previous:
        st.session_state["contador"] -= 1
    
    if st.session_state["contador"] < 0:
        st.session_state["contador"] = st.session_state["totalgraficos"] + st.session_state["contador"]
    if st.session_state["contador"] > (st.session_state["totalgraficos"]-1):
        st.session_state["contador"] = st.session_state["contador"]-st.session_state["totalgraficos"]

 

    with ColB:
       
        # Contador = 1 --> representacion hla
        # Contador = 2 --> representacion tumores
        # Contador = 3 --> Anchor non anchor pos neg
        # Contador = 4 --> Cantidad de mutaciones por posicion
        # Contador = 5 --> Matriz frecuencia long 9
        # Contador = 6 --> Matriz frecuencia long 10
        # contador = 7 --> cant long 9 y long 10 y dividir pos y neg
        container = st.container(border=True)
        with container:
    
            if (st.session_state["contador"]+1) == 1:   
                st.subheader("HLA Alleles")         
                st.plotly_chart(fig_HLA, theme="streamlit", use_container_width=True)
                st.caption('Neoantigen frequency by HLA allele, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each HLA allele.')
            if (st.session_state["contador"]+1) == 2:
                st.subheader("Tumor types")
                st.plotly_chart(fig_tumor, theme="streamlit", use_container_width=True)
                st.caption('Pie chart showing the distribution of tumor types in the database.')
            if (st.session_state["contador"]+1) == 3:
                st.subheader("Mutation types")
                st.plotly_chart(fig_PositionType_NeoType, use_container_width=True)
                st.caption('Bar chart showing the frequency distribution of anchor and non-anchor position types, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each position type.')
            if (st.session_state["contador"]+1) == 4:
                st.subheader("Mutation's position")
                st.plotly_chart(fig_mutPosition_NeoType, use_container_width=True)
                st.caption('Bar chart showing the frequency of mutations at each position, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each position.')
            if (st.session_state["contador"]+1) == 5:
                st.subheader("Neoantigens' length")
                st.plotly_chart(fig_Length_NeoType, use_container_width=True)
                st.caption('Bar chart showing the frequency of neoantigen lengths, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each length.')
            if (st.session_state["contador"]+1) == 6:
                st.subheader("Length 9 neoantigens' PFM matrices")
                st.plotly_chart(fig_mat_frec_9, use_container_width=True)
                st.caption('PFM matrices of 9-mer neoantigens: A. Immunogenic, B. Non-immunogenic, C. Difference. Each matrix shows the frequency of each amino acid at each position.')
            if (st.session_state["contador"]+1) == 7:   
                st.subheader("Length 10 neoantigens' PFM matrices")
                st.plotly_chart(fig_mat_frec_10, use_container_width=True)
                st.caption('PFM matrices of 10-mer neoantigens: A. Immunogenic, B. Non-immunogenic, C. Difference. Each matrix shows the frequency of each amino acid at each position.')





with tab2:

    #Creación de listas con valores unicos y ordenadas alfabeticamente    
    Tumor_list = df_ITSNdb_complete.Tumor.unique()  
    Tumor_list = Tumor_list.tolist()
    Tumor_list.sort()
    
    Neoag_list = df_ITSNdb_complete.Neoantigen.unique() 
    Neoag_list = Neoag_list.tolist()
    Neoag_list.sort()

    wt_list = df_ITSNdb_complete.WT.unique()
    wt_list = wt_list.tolist()
    wt_list.sort()

    HLA_gene_list = df_ITSNdb_complete.HLA_GENE.unique()
    HLA_gene_list = HLA_gene_list.tolist()
    HLA_gene_list.sort()

    HLA_Allele_list = df_ITSNdb_complete.HLA_ALLELE_GROUP.unique()
    HLA_Allele_list = HLA_Allele_list.tolist()
    HLA_Allele_list.sort()

    HLA_protein_list = df_ITSNdb_complete.HLA_PROTEIN.unique()
    HLA_protein_list = HLA_protein_list.tolist()
    HLA_protein_list.sort()

    Gene_list = df_ITSNdb_complete.GeneSymbol.unique()
    Gene_list = Gene_list.tolist()
    Gene_list.sort()

    Length_list = df_ITSNdb_complete.Length.unique()
    Length_list = Length_list.tolist()
    Length_list.sort()
    Length_list = ["Any"] + Length_list

    Mutpos_list = df_ITSNdb_complete.mutPosition.unique()
    Mutpos_list = Mutpos_list.tolist()
    Mutpos_list.sort()


    #Creacion de objetos interctivos 
    # label_expander = "**Filters**"
    # with st.expander(label = label_expander,
    #                   expanded = True):
        
    
    # C1, C2, C3, C4 = st.columns(4)

    # with C1:
    #     Neoag_multiselect = st.multiselect(label = "Neoantigens", 
    #                                 options = Neoag_list)
    #     gene_multiselect = st.multiselect (label = "Gene Symbol",
    #                                 options = Gene_list)
    #     Neotype_list=st.radio(label = "Inmunogenecity",
    #                             options = ["Both","Positive", "Negative"], 
    #                             index=0)

    # with C2:
    #     Wt_multiselect = st.multiselect(label = "WildType",
    #                                     options = wt_list)
    #     HLA_gene_multiselect = st.multiselect(label = "HLA gene",
    #                                     options = HLA_gene_list)
    #     Pos_type=st.radio(label = "Position type",
    #                     options = ["Both","Anchor", "Non-anchor"],
    #                     index=0)
    
    # with C3:
    #     mutpos_multiselect = st.multiselect(label = "Mutation position",
    #                                 options = Mutpos_list)
    #     HLA_alllele_multiselect = st.multiselect(label = "HLA allele group",
    #                                         options = HLA_Allele_list)
    #     Length_radio = st.radio(label = "Length", 
    #                     options = Length_list, 
    #                     index = 0)


    # with C4:
    #     Tumor_multiselect = st.multiselect(label = "Tumors", 
    #                                     options = Tumor_list) 
    #     HLA_protein_multiselect = st.multiselect(label = "HLA protein",
    #                                         options = HLA_protein_list)
    #     # hla_nom = st.popover(label="HLA nomenclature", use_container_width=True )
    #     # with hla_nom:
    #     #     st.image("HLA_NOM.png") 
    #     st.image("HLA_NOM.png")

    with st.expander(label = '**Filters:** ',
                      expanded = True):
        C1, C2, C3, C4= st.columns(4)

        with C1:
            Neoag_multiselect = st.multiselect(label = "Neoantigens", 
                                        options = Neoag_list)
        with C2:
            Wt_multiselect = st.multiselect(label = "WildType",
                                            options = wt_list)
        with C3:
            mutpos_multiselect = st.multiselect(label = "Mutation position",
                                        options = Mutpos_list)
            
        with C1:
            HLA_gene_multiselect = st.multiselect(label = "HLA gene",
                                            options = HLA_gene_list)
        with C2:
            HLA_alllele_multiselect = st.multiselect(label = "HLA allele group",
                                                options = HLA_Allele_list)
        with C3:
            HLA_protein_multiselect = st.multiselect(label = "HLA protein",
                                                options = HLA_protein_list)
        with C4:
            st.image("HLA_NOM.png", width= 300)    

        c13, c23, c33 = st.columns([1,1,2])

        with c13:
            gene_multiselect = st.multiselect (label = "Gene Symbol",
                                        options = Gene_list)
        with c23:
            Tumor_multiselect = st.multiselect(label = "Tumors", 
                                            options = Tumor_list)
        with c33:

            C12, C22, C32, C42 = st.columns(4)

            with C12:
                Neotype_list=st.radio(label = "Inmunogenecity",
                        options = ["Both","Positive", "Negative"], 
                        index=0)
                        

            with C22:
                Pos_type=st.radio(label = "Position type",
                        options = ["Both","Anchor", "Non-anchor"],
                        index=0)
            with C32:
                Length_radio = st.radio(label = "Length", 
                                options = Length_list, 
                                index = 0)


       




            




    #Ajuste para elementos vacios

    if not gene_multiselect:
        sel_gene = Gene_list
    else: 
        sel_gene = gene_multiselect

    if not mutpos_multiselect:
        mutpos_sel = Mutpos_list
    else:
        mutpos_sel = mutpos_multiselect

    if str(Length_radio) == "Any":
        Length_sel = Length_list
    else:
        Length_sel = [Length_radio]

    if Neotype_list=="Positive":
        Neotype_sel=["Positive"]
    elif Neotype_list=="Negative":
        Neotype_sel=["Negative"]
    else:
        Neotype_sel=["Positive","Negative"]
    

    if not HLA_protein_multiselect: 
        sel_hla_protein = HLA_protein_list
    else:
        sel_hla_protein = HLA_protein_multiselect

    if not HLA_alllele_multiselect:
        sel_hla_allele = HLA_Allele_list
    else:
        sel_hla_allele = HLA_alllele_multiselect

    if not HLA_gene_multiselect:
        sel_hla_gene = HLA_gene_list
    else: 
        sel_hla_gene = HLA_gene_multiselect

    if not Wt_multiselect:
        sel_wt = wt_list
    else:
        sel_wt = Wt_multiselect

    if not Neoag_multiselect:
        sel_neoag = Neoag_list
    else:
        sel_neoag = Neoag_multiselect

    if not Tumor_multiselect:
        sel_tumor = Tumor_list
    else:
        sel_tumor = Tumor_multiselect


    
    if Pos_type=="Anchor":
        Pos_type_sel=["Anchor"]
    elif Pos_type=="Non-anchor":
        Pos_type_sel=["Non-anchor"]
    else:
        Pos_type_sel=["Anchor","Non-anchor"]    

    #Filtrado 

    df_filtered = df_ITSNdb_complete[df_ITSNdb_complete["Tumor"].isin(sel_tumor)
                                     & df_ITSNdb_complete["Neoantigen"].isin(sel_neoag)
                                     & df_ITSNdb_complete["WT"].isin(sel_wt)
                                     & df_ITSNdb_complete["HLA_GENE"].isin(sel_hla_gene)
                                     & df_ITSNdb_complete["HLA_ALLELE_GROUP"].isin(sel_hla_allele)
                                     & df_ITSNdb_complete["HLA_PROTEIN"].isin(sel_hla_protein)
                                     & df_ITSNdb_complete["NeoType"].isin(Neotype_sel)
                                     & df_ITSNdb_complete["GeneSymbol"].isin(sel_gene)
                                     & df_ITSNdb_complete["Length"].isin(Length_sel)
                                     & df_ITSNdb_complete["mutPosition"].isin(mutpos_sel)
                                     & df_ITSNdb_complete["PositionType"].isin(Pos_type_sel)]

    #st.write("   ")
    
    df_filtered["Select"] = 0
    
    # with st.expander(label="Database", expanded=True):
    st.subheader("ITSN database")
    df_filtered = st.data_editor(data= df_filtered,
                                column_config = {
                                    "Select": st.column_config.CheckboxColumn("Select",
                                                                            help="Select the files you want to see in the 'Neoantigen information' section", 
                                                                            default=False),
                                    "Author": None,
                                    "Paper": None,
                                    "HLA_GENE": None,
                                    "HLA_ALLELE_GROUP": None, 
                                    "HLA_PROTEIN": None,
                                },
                                column_order = ("Select", "Neoantigen", "WT", "NeoType", "mutPosition", "PositionType",  "Length", "HLA", "GeneSymbol", "Tumor"),
                                disabled=["Neoantigen", "WT", "NeoType", "mutPosition", "Length", "HLA", "GeneSymbol", "Tumor"],
                                hide_index=True,
                                use_container_width = True,
                                height = 400
                                    ) 
    

#with tab3:

    st.subheader("Neoantigen information")

    df_Selected = df_filtered[df_filtered["Select"]==1]



    if df_Selected.empty:
        st.caption("Please, select neoantigens in the database above to see more detailed information about particular neoantigens")
    #st.dataframe(df_Selected)
    for i in df_Selected.index:
            
        #EXTRACCIÓN DE LOS DATOS DEL REGISTRO SELECCIONADO            
        neo_s=df_Selected.loc[i,"Neoantigen"]
        wt_s=df_Selected.loc[i,"WT"]
        author_s=df_Selected.loc[i,"Author"]
        paper_s=df_Selected.loc[i,"Paper"]
        gene_s=df_Selected.loc[i,"GeneSymbol"]
        length_s= df_Selected.loc[i, "Length"]
        mut_pos_s = df_Selected.loc[i, "mutPosition"]
        #Genera lista a partir de string
        list_Neo = list(neo_s)
        list_Wt = list(wt_s)
        list_espacios = []
        lista_columnas = []
        for i in range(0, len(list_Neo)):
            if list_Neo[i] == list_Wt[i]:
                list_espacios.append("🟢")
            else:
                list_espacios.append("❌")
            lista_columnas.append(i+1)       

        #Lista de genes en la base de datos
        lista_genes = df_ITSNdb_complete['GeneSymbol'].tolist()
        lista_neoantigenos = df_ITSNdb_complete['Neoantigen'].tolist()
        genes_muy_pesados = ['KMT2C', 'DYNC1H1', 'HUWE1', 'ALMS1', 'UTRN']
        genes_predichos_txt = 'Genes_predichos.txt'
        with open(genes_predichos_txt, 'r') as f:
            lista_wt_predichos = [linea.strip() for linea in f]
        lista_sin_predicción = []
        #Ubicar predicciones de wt faltantes
        for gen in lista_genes:
            gen_nuevo = gen.lower()+'_wt'
            if gen_nuevo not in lista_wt_predichos:
                lista_sin_predicción.append(gen)
        
        #generar rutas a los archivos.lower()+'_'
        path_pdb_wt_git = 'https://raw.githubusercontent.com/roflesia/ITSN_App/main/PDB/WT/'
        path_pdb_mut_git = 'https://raw.githubusercontent.com/roflesia/ITSN_App/main/PDB/MUT/'
        
        path_pae_wt_git = 'https://bitbucket.org/itsn_app/itsn_app/raw/main/PAE_WT/'
        path_pae_mut_git = 'https://bitbucket.org/itsn_app/itsn_app/raw/main/PAE_MUT/'

        path_fasta_wt_git = 'https://raw.githubusercontent.com/roflesia/ITSN_App/main/WT_FASTA/'
        path_fasta_mut_git = 'https://raw.githubusercontent.com/roflesia/ITSN_App/main/NEO_FASTAS/'

        gen_minu = gene_s.lower()
        neoag_minu = neo_s.lower()         

        #SE GENERA UN EXPANDER PARA CADA REGISTRO SELECCIONADO
        with st.expander(gene_s+" - "+neo_s, expanded=True):
            #URL gencard
            url_gc_ic = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
            url_up_ic = 'https://www.uniprot.org/uniprotkb/'
            #gencard name
            genecard_id = ''
            lista_genes = df_ITSNdb_GENES['GeneSymbol'].tolist()
            lista_genescard =df_ITSNdb_GENES['Gene_Card'].tolist()
            lista_unitprot=df_ITSNdb_GENES['Uniprot_code'].tolist()
            lista_seq_len=df_ITSNdb_GENES['Seq_Lenght'].tolist()
            for gen in range(0, len(lista_genes)):
                if gene_s == lista_genes[gen]:
                    genecard_id = lista_genescard[gen]
                    unitprot_id = lista_unitprot[gen]
                    seq_len =lista_seq_len[gen]

            ptm_df = df_ITSNdb_GENES[df_ITSNdb_GENES['Neoantigen'] == neo_s]

            if not ptm_df.empty:
                # Obtener los valores de PTMscore para WT y MUT
                ptm_wt = ptm_df['PTMscore_wt'].values[0]
                ptm_mut = ptm_df['PTMscore_mut'].values[0]
                seq_mut = ptm_df['GENE_SEQ_mut'].values[0]

            else:
                ptm_wt = None
                ptm_mut = None
                st.write("Error prediccion ptm")


            secuencia_completa = seq_mut
            neoantigeno = neo_s
            posicion_mutada = mut_pos_s  # Posición 1-indexed dentro del neoantígeno

            posicion_en_secuencia_completa = encontrar_posicion_mutacion(secuencia_completa, neoantigeno, posicion_mutada)

            
            listas = [list_Neo, list_espacios, list_Wt]
            df_compared_mut = pd.DataFrame(columns= lista_columnas, index=["NeoAntigen", "", "WildType"] , data=listas)
            
            t1,t2, t3= st.columns([2,3,2])
            with t1:
                st.metric(label = ("Neoantigen:"),
                        value = neo_s)
            
                st.metric(label = ("Wildtype:"),
                        value = wt_s)
                
            with t2:
                r1, r2 = st.columns(2)
                with r1:


                    st.write("	:straight_ruler:", "**Protein lenght:** ", int(seq_len))
                    st.write("	:straight_ruler:", "**Neoantigen lenght:** ", length_s)
                with r2:
                    if posicion_en_secuencia_completa: 
                        st.write(":pushpin:", "**Mutation position on protein:** ", posicion_en_secuencia_completa)
                    else:
                        st.write(":pushpin:", "**Mutation position on protein: X** ")
                    st.write(":pushpin:", "**Mutation position:** ", mut_pos_s)
                st.dataframe(df_compared_mut)
                

            with t3:
                # #URL gencard
                # url_gc_ic = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
                # url_up_ic = 'https://www.uniprot.org/uniprotkb/'
                # #gencard name
                # genecard_id = ''
                # lista_genes = df_ITSNdb_GENES['GeneSymbol'].tolist()
                # lista_genescard =df_ITSNdb_GENES['Gene_Card'].tolist()
                # lista_unitprot=df_ITSNdb_GENES['Uniprot_code'].tolist()
                # lista_seq_len=df_ITSNdb_GENES['Seq_Lenght'].tolist()
                
                # for gen in range(0, len(lista_genes)):
                #     if gene_s == lista_genes[gen]:
                #         genecard_id = lista_genescard[gen]
                #         unitprot_id = lista_unitprot[gen]
                #         seq_len =lista_seq_len[gen]
                
                url_gc = url_gc_ic + genecard_id.strip()
                url_up = url_up_ic + unitprot_id.strip()                
                label1 = ('GeneCard: ' + genecard_id)
                label2 = ('UniProt: '+ unitprot_id)


                st.write(" **Gene:** ", gene_s)

                st.page_link(page = url_gc, label = label1, icon="🧬")              
                st.page_link(page = url_up, label = label2, icon="🔬")                        
                st.write("**Link to  the author´s reference paper:** ")
                st.page_link(page = paper_s, label = author_s, icon="📑")

            #ptm Score
            # Filtrar el DataFrame por el neoantígeno proporcionado
            # ptm_df = df_ITSNdb_GENES[df_ITSNdb_GENES['Neoantigen'] == neo_s]

            # if not ptm_df.empty:
            #     # Obtener los valores de PTMscore para WT y MUT
            #     ptm_wt = ptm_df['PTMscore_wt'].values[0]
            #     ptm_mut = ptm_df['PTMscore_mut'].values[0]
            #     seq_mut = ptm_df['GENE_SEQ_mut'].values[0]

            # else:
            #     ptm_wt = None
            #     ptm_mut = None
            #     st.write("Error prediccion ptm")

            # st.write("PTM Score WT: ", ptm_wt)
            # st.write("PTM Score MUT: ", ptm_mut)

            # Ejemplo de uso:
            
            # secuencia_completa = seq_mut
            # neoantigeno = neo_s
            # posicion_mutada = mut_pos_s  # Posición 1-indexed dentro del neoantígeno

            # posicion_en_secuencia_completa = encontrar_posicion_mutacion(secuencia_completa, neoantigeno, posicion_mutada)
            # if posicion_en_secuencia_completa:
            #     st.write("La mutación está en la posición ", posicion_en_secuencia_completa, " de la secuencia completa.")
            # else:
            #     st.write("El neoantígeno no se encuentra en la secuencia completa.")
        # #EXTRACCIÓN DE LOS DATOS DEL REGISTRO SELECCIONADO            
        
        
            if gene_s not in lista_sin_predicción:
                
                with st.form(gene_s+neo_s):
                    st.subheader("AlphaFold 3 prediction")
                    Sel_pred =st.radio(label= "Pred", options=["Wildtype vs. Mutation comparative", "Mutation results", "Wildtype results"], index=0, key=gene_s+neo_s+"2", horizontal=True, captions=None, label_visibility="collapsed")
                    sel_superimpose = st.radio(
                            label = "Select superimpose mode",
                            options = [ "PLDDT > 70 ", "Global"], 
                            key = neoag_minu,
                            horizontal=True)                              
                    # with hilera2:
                    #     st3d_view = st.checkbox(label= "3D Structure", value=False, key= "3D Structure"+neo_s, label_visibility="visible")
                    #     pae_view = st.checkbox(label= "PAE", value=False, key= "pae"+neo_s, label_visibility="visible")
                    #     plddt_view = st.checkbox(label= "Confidence Scores", value=False, key= "PLDDT"+neo_s, label_visibility="visible")
                    

                    # Every form must have a submit button.
                    submitted = st.form_submit_button("Submit")

            if submitted:
                if seq_len > 3000:
                    st.info('Sequences longer than 2000 AA are too big to upload some graphs', icon="ℹ️")
                elif seq_len > 1000:
                    st.info('Larger sequences might take a moment to load', icon="ℹ️")
                #Obtener PDB
                path_pdb_wt = path_pdb_wt_git + gen_minu + '_wt.pdb' 
                path_pdb_mut = path_pdb_mut_git + gen_minu + '_' + neoag_minu + '.pdb'
                path_fasta_wt = path_fasta_wt_git + gene_s + '_wt.fasta' 
                path_fasta_mut = path_fasta_mut_git + gene_s + '_' + neo_s + '.fasta'
                df_pdb_wt = GET_PDB_CA(path_pdb_wt)
                df_pdb_mut = GET_PDB_CA(path_pdb_mut)

                with urllib.request.urlopen(path_pdb_wt) as response:
                        wt_pdb = response.read().decode('utf-8')
                with urllib.request.urlopen(path_pdb_mut) as response:
                        mut_pdb = response.read().decode('utf-8')



                l1, l2 = st.columns([3,2])
                with l1:
                    st.header('3D structure')
                      
                    if sel_superimpose == "Global":
                        sup_pdb =  superimpose_proteins(wt_pdb, mut_pdb)
                    else:
                        sup_pdb =  superimpose_pdb_with_threshold(path_pdb_wt, path_pdb_mut)

                    Spin  = False
                    
                #Obtener PAE
                if  gene_s not in genes_muy_pesados:       
                    path_pae_wt = path_pae_wt_git + gen_minu + '_wt.json' 
                    path_pae_mut = path_pae_mut_git + gen_minu + '_' + neoag_minu + '.json'
                    df_pae_wt = GET_PAE_DF(path_pae_wt)
                    df_pae_mut = GET_PAE_DF(path_pae_mut)
                
                view_3d = py3Dmol.view(width=800, height=500) 


            
                
                if Sel_pred == "Wildtype vs. Mutation comparative":
                    st.write("Wildtype vs. Mutation comparative")
                    
                    if  gene_s not in genes_muy_pesados:
                        fig_pae = GET_PAE_DIF(df_pae_wt, df_pae_mut)
                    
                    view_3d.addModel(wt_pdb, 'pdb')
                    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0x51adbe'}})  # Modelo 0 con color azul
                    view_3d.addModel(sup_pdb, 'pdb')
                    view_3d.setStyle({'model': 1}, {"cartoon": {'color': '0xdb7093'}})  # Modelo 0 con color azul

                elif Sel_pred == "Mutation results":
                    st.write("Mutation results")
                    if  gene_s not in genes_muy_pesados:
                        fig_pae = GET_PAE_GRAPH(df_pae_mut, posicion_en_secuencia_completa)
                    view_3d.addModel(mut_pdb, 'pdb')
                    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0xdb7093'}})  # Modelo 0 con color azul                            
                else:
                    st.write("Wildtype results")
                    if  gene_s not in genes_muy_pesados:
                        fig_pae = GET_PAE_GRAPH(df_pae_wt, posicion_en_secuencia_completa)
                    view_3d.addModel(wt_pdb, 'pdb')
                    view_3d.setStyle({'model': 0}, {"cartoon": {'color': '0x51adbe'}})  # Modelo 0 con color azul
                
                
                with l1:
                    
                    view_3d.spin(Spin)
                    view_3d.zoomTo()
                    #showmol(view_3d)
                    st.components.v1.html(view_3d._make_html(), height=500, width=800)

                with l2:
                    st.header('PAE')
                    if gene_s not in genes_muy_pesados:
                        st.plotly_chart(fig_pae)
                    else:
                        st.warning("This structure is too large to graph its PAE.")

            
                ## aa distance
                
                reference_pdb_path = wt_pdb
                sample_pdb_path = sup_pdb
                distancias, tmscore, rmsd = aadistance(reference_pdb_path, sample_pdb_path, posicion_en_secuencia_completa)
                ##PLDDT
                ref_wt = "WildType"
                ref_mut = "Mutation"
                fig_plddt = GET_PLDDTS(df_pdb_wt, ref_wt, df_pdb_mut, ref_mut, posicion_en_secuencia_completa, Sel_pred)


                if Sel_pred == "Wildtype vs. Mutation comparative":
                    w1, w2, w3 = st.columns([1,1,3])
                    with w1:
                        st.metric(label = ("PTM Score WT:"),
                            value = ptm_wt)
                    with w2:
                        st.metric(label = ("PTM Score MUT:"),
                            value = ptm_mut)
                elif Sel_pred == "Mutation results":
                    st.metric(label = ("PTM Score MUT:"),
                            value = ptm_mut)
                else:
                    st.metric(label = ("PTM Score WT:"),
                             value = ptm_wt)
                
                #Grafico plddt
                st.header('PLDDT')
                st.plotly_chart(fig_plddt)


                if Sel_pred == "Wildtype vs. Mutation comparative":
                    st.header('AA distance')
                    u1, u2, u3 = st.columns([1,1,3])
                    with u1:
                        st.metric(label = ("TM Score:"),
                            value = round(tmscore, 3))
                    with u2:
                        st.metric(label = ("RMSD:"),
                            value = (str(round(rmsd,3)) + " Å"))
                    st.plotly_chart(distancias)
                

                
                

                
                json_url = "https://github.com/tu_usuario/tu_repo/ruta_del_archivo.json"
                fasta_url = "https://bitbucket.org/tu_usuario/tu_repo/ruta_del_archivo.fasta"

                # path_pdb_wt 
                # path_pdb_mut 
                # path_pae_wt
                # path_pae_mut
                # path_fasta_wt
                # path_fasta_mut
                def descargar_archivo(url):
                    response = requests.get(url)
                    response.raise_for_status()  # Asegura que no haya errores en la descarga
                    return response.content

                dl1, dl2, dl3, dl4, dl5, dl6 = st.columns(6)
                with dl1:
                    # Botón de descarga para archivo .pdb
                    pdb_content = descargar_archivo(path_pdb_wt)
                    st.download_button(label="Download wildtype's pdb", data=pdb_content, file_name= gene_s + "_wt.pdb", key = "pdb" + path_pdb_wt)
                with dl2:
                    pdb_content = descargar_archivo(path_pdb_mut)
                    st.download_button(label="Download mutation's pdb", data=pdb_content, file_name= gene_s + '_' + neo_s + ".pdb", key = "pdb" + path_pdb_mut)
                with dl3:
                    # Botón de descarga para archivo .json
                    json_content = descargar_archivo(path_pae_wt)
                    st.download_button(label="Download wildtype's pae", data=json_content, file_name= gene_s + "_wt.json", key = "pae" + path_pae_wt)
                with dl4:
                    json_content = descargar_archivo(path_pae_mut)
                    st.download_button(label="Download mutation's pae", data=json_content, file_name= gene_s + '_' + neo_s + ".json", key= "pae" + path_pae_mut)
                # Botón de descarga para archivo .fasta
                with dl5:
                    fasta_content = descargar_archivo(path_fasta_wt)
                    st.download_button(label="Download wildtype's fasta", data=fasta_content, file_name= gene_s + "_wt.fasta", key= "fasta" + path_fasta_wt)
                with dl6:
                    fasta_content = descargar_archivo(path_fasta_mut)
                    st.download_button(label="Download mutation's fasta", data=fasta_content, file_name= gene_s + '_' + neo_s + ".fasta", key = "fasta" + path_fasta_mut)


with tab3:    
    ITSNdb = pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/ITSNdb.csv", index_col=False)
    Val_ds = pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/Val_dataset.csv", index_col=False)
    TNB_ds = pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/TNB_dataset.csv", index_col=False)

    #Lista de columnas 

    list_col_itsndb_or = ITSNdb.columns.to_list()
    list_col_valds_or = Val_ds.columns.to_list()
    list_col_tnb_or  = TNB_ds.columns.to_list()

    list_col_itsndb_or.append("Prediction_result")
    list_col_valds_or.append("Prediction_result")
    list_col_tnb_or.append("Prediction_result")


    @st.cache_data
    def convert_df(df):
        # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode("utf-8")

    ITSNdb_csv = convert_df(ITSNdb)
    Val_ds_csv = convert_df(Val_ds)
    TNB_ds_csv = convert_df(TNB_ds)

    ITSNdb_pred_upload = pd.DataFrame()
    Valds_pred_upload = pd.DataFrame()
    TNBds_pred_upload = pd.DataFrame()

    #Funcion roc curve metrics 
    def get_ROC_metrics(thr_array, fpr_array, tpr_array, pos_count, neg_count):
        Sp_array, Se_array, ppv_array, npv_array, tp_array, fp_array, tn_array, fn_array = [], [], [], [], [], [], [], []
        dop_array = []
        F1_Score_array = []
        for i in range(0, len(thr_array)):
            Sp_array.append(1 - fpr_array[i])
            Se_array.append(tpr_array[i])
            
            tp_array.append(tpr_array[i] * pos_count)
            tn_array.append((1-fpr_array[i]) * neg_count)
            fp_array.append(fpr_array[i] * neg_count)
            fn_array.append((1-tpr_array[i]) * pos_count)

            if (tpr_array[i] * pos_count + fpr_array[i] * neg_count) == 0:
                ppv_array.append(None)
            else:
                ppv_array.append(tpr_array[i] * pos_count / (tpr_array[i] * pos_count + fpr_array[i] * neg_count))
            if (1 - fpr_array[i]) * neg_count + (1 - tpr_array[i]) * pos_count == 0:
                npv_array.append(None)
            else:
                npv_array.append((1 - fpr_array[i]) * neg_count / ((1 - fpr_array[i]) * neg_count + (1 - tpr_array[i]) * pos_count))

            if ppv_array[i] == None or npv_array[i]  == None:
                dop_array.append(1500)
            else:
                dop_array.append(np.sqrt((Se_array[i] - 1) ** 2 + (Sp_array[i] - 1) ** 2 + 
                                        (ppv_array[i] - 1) ** 2 + (npv_array[i] - 1) ** 2))

            f1score_value = tp_array[i]/(tp_array[i]+1/2*(fp_array[i]+fn_array[i]))
            F1_Score_array.append(f1score_value)



        min_index = 1
        for j in range(1, len(F1_Score_array)):
            if dop_array[min_index] > dop_array[j]:
                    min_index = j


        dict_1 = {
            "Threshold" : thr_array,
            "Specificity" : Sp_array,
            "Sensitivity" : Se_array,
            "TPR" : tpr_array,
            "FPR" : fpr_array,
            "TP": tp_array,
            "TN" : tn_array,
            "FP" : fp_array,
            "FN" : fn_array,
            "PPV" : ppv_array,
            "NPV": npv_array,
            "DOP": dop_array,
            "F1 Score": F1_Score_array
        }
        df_roc_metrics = pd.DataFrame(dict_1)

        min_dop = (dop_array[min_index])
        dop_thr = (thr_array[min_index])
        dop_f1 = (F1_Score_array[min_index])

        auc = metrics.auc(new_fpr, new_tpr)
        

        return df_roc_metrics, min_dop, dop_thr, dop_f1, auc
  #st.header("Software Validation")

    with st.expander("", expanded=True):

        c51, c52 = st.columns([3,5])

        with c51:
            st.image("Procedimineto_carga.png")
            st.caption("*Create a new column in each file to store the predictions. The name of this column should be “Prediction_result” and it must be the last column of each file.")
            st.caption("**The names of the csv files must be 'ITSNdb_pred.csv', 'Valds_pred.csv' and 'TNBds_pred.csv'")
            
        with c52:
            st.subheader("1. Download the datasets")
            cb1, cb2, cb3 = st.columns(3)
            with cb1:
                st.download_button(
                    label="Download ITSNdb",
                    data=ITSNdb_csv,
                    file_name="ITSNdb.csv",
                    mime="text/csv", 
                    use_container_width= True)
            with cb2:
                st.download_button(
                    label="Download Validation DataSet",
                    data=Val_ds_csv,
                    file_name="Valds.csv",
                    mime="text/csv", 
                    use_container_width= True)
            with cb3:
                st.download_button(
                    label="Download TNB DataSet",
                    data=TNB_ds_csv,
                    file_name="TNBds.csv",
                    mime="text/csv", 
                    use_container_width= True)
            
            st.subheader("2. Make the predictions")

            st.subheader("3. Upload the datasets with the prediction results")
            uploaded_files = st.file_uploader("Choose a CSV file", accept_multiple_files=True)
            flag_ITSNdb = 0
            flag_Valds = 0
            flag_TNB = 0

            flag_col_ITSNdb = 0
            flag_col_valds = 0
            flag_col_tnb = 0
            for uploaded_file in uploaded_files:
                if uploaded_file.name == "ITSNdb_pred.csv":
                    flag_ITSNdb = 1
                    ITSNdb_pred_upload = pd.read_csv(uploaded_file)
                    list_col_itsndb = ITSNdb_pred_upload.columns.to_list()
                    if set(list_col_itsndb_or) == set(list_col_itsndb):
                        flag_col_ITSNdb = 1
                if uploaded_file.name == "Valds_pred.csv":
                    Valds_pred_upload = pd.read_csv(uploaded_file)
                    flag_Valds = 1
                    list_col_valds = Valds_pred_upload.columns.to_list()
                    if set(list_col_valds) == set(list_col_valds_or):
                        flag_col_valds = 1
                if uploaded_file.name == "TNBds_pred.csv":
                    flag_TNB = 1
                    TNBds_pred_upload = pd.read_csv(uploaded_file)     
                    list_col_tnb = TNBds_pred_upload.columns.to_list()
                    if set(list_col_tnb) == set(list_col_tnb_or):
                        flag_col_tnb = 1     

            st.subheader("Based on the software's predicted values, is there a positive correlation between these values and immunogenicity?")
            correlation = st.radio(label= "Based on the software's predicted values, is there a positive correlation between these values and immunogenicity?", label_visibility= "collapsed", options = ["Yes", "No"], captions= ["Neoantigens with values above the threshold are considered immunogenic.", "Neoantigens with values above the threshold are considered non-immunogenic."])

    flag_ag = 0
    flag_ag = flag_ITSNdb and flag_Valds and flag_TNB and flag_col_ITSNdb and flag_col_valds and flag_col_tnb


    if not flag_ag:
        if not (flag_ITSNdb or flag_Valds or flag_TNB):
            st.warning("No files uploaded")
        else: 
            if not (flag_ITSNdb): 
                st.warning("No ITSNdb_pred.csv file uploaded")
            elif not (flag_col_ITSNdb):
                st.error("ITSNdb_pred.csv´s columns names do not match the expected names")
            if not (flag_Valds): 
                st.warning("No Valds_pred.csv file uploaded") 
            elif not (flag_col_valds):
                st.error("Valds_pred.csv´s columns names do not match the expected names")
            if not (flag_TNB): 
                st.warning("No TNBds_pred.csv file uploaded")          
            elif not (flag_col_tnb):
                st.error("TNBds_pred.csv´s columns names do not match the expected names")
    else:   
        st.success("All files were succefully uploaded")

        #st.write(correlation)
        
        st.header("ITSNdb")

        Pos_count_ITSNdb =  ITSNdb.value_counts(ITSNdb['NeoType']).iat[0]
        Neg_Count_ITSNdb = ITSNdb.value_counts(ITSNdb['NeoType']).iat[1]
            
        ITSNdb_pred_upload["NeoType_01"] = 0
        for i in range(len(ITSNdb_pred_upload)):
            if ITSNdb_pred_upload['NeoType'].iloc[i] == "Positive":
                ITSNdb_pred_upload["NeoType_01"].iloc[i] = 1
            if ITSNdb_pred_upload['NeoType'].iloc[i] == "Negative":
                ITSNdb_pred_upload["NeoType_01"].iloc[i] = 0
        
        #Curva roc
        #aplica funcion 
        
        pos_lb= 0
        if correlation == "Yes":
            pos_lb = 0
        else:
            pos_lb = 1



        roc = ru.compute_roc(X=ITSNdb_pred_upload["Prediction_result"], y=ITSNdb_pred_upload["NeoType_01"], pos_label=pos_lb )
        #Valores intermedios
        thr_roc = roc.thr
        fpr_roc = roc.fpr 
        tpr_roc = roc.tpr
        auc_roc = roc.auc

        new_tpr = []
        new_fpr = []
        for i in range(len(thr_roc)):
            new_fpr.append(1-fpr_roc[i])
            new_tpr.append(1-tpr_roc[i])

        

        st.write("Roc curve metrics")
        df_roc_metrics, min_dop, dop_thr, dop_f1, auc_roc = get_ROC_metrics(thr_array= thr_roc, fpr_array= new_fpr, tpr_array= new_tpr, pos_count= Pos_count_ITSNdb, neg_count= Neg_Count_ITSNdb)
        #st.dataframe(df_roc_metrics)

        d1, d2 = st.columns([1,4])
        with d1:
            st.metric(label = ("Min DOP:"),
                            value = round(min_dop, 4))
        
            st.metric(label = ("DOP THR:"),
                            value = round(dop_thr, 4))
        
            st.metric(label = ("DOP F1:"),
                            value =round(dop_f1, 4))
    
            st.metric(label = ("AUC:"),
                            value = round(auc_roc, 4))    
        style_metric_cards()


        Roc_Curve = px.area(
            x=new_fpr, y=new_tpr,
            labels=dict(x='1-Specificity', y='Sensitivity'),
            width=500, height=500
        )
        Roc_Curve.add_shape(
            type='line', line=dict(dash='dash'),
            x0=0, x1=1, y0=0, y1=1
        )

        Roc_Curve.update_yaxes(scaleanchor="x", scaleratio=1)
        Roc_Curve.update_xaxes(constrain='domain')

        with d2:    
            st.plotly_chart(Roc_Curve, use_container_width=True )

        with st.expander(label = "ROC Metrics dataframe", expanded = False):
            st.dataframe(df_roc_metrics, use_container_width=True)        


        st.header("Validation DataSet")

        #df_val_pred = pd.read_csv("Val_dataset_ent_MHC_results.csv")


        #crear una columna nueva con la prediccion segun el threshold

        Valds_pred_upload["pred_result"] = 0

        if pos_lb == 0:
            for i in range(len(Valds_pred_upload)):
                if Valds_pred_upload["Prediction_result"].iloc[i] < dop_thr:
                    Valds_pred_upload["pred_result"].iloc[i] = "Negative"
                else: 
                    Valds_pred_upload["pred_result"].iloc[i] = "Positive"
        if pos_lb == 1:
            for i in range(len(Valds_pred_upload)):
                if Valds_pred_upload["Prediction_result"].iloc[i] >= dop_thr:
                    Valds_pred_upload["pred_result"].iloc[i] = "Negative"
                else: 
                    Valds_pred_upload["pred_result"].iloc[i] = "Positive"


        #contamos tp fp tn y fn 
        valds_tp = 0
        valds_fp = 0
        valds_tn = 0
        valds_fn = 0

        Valds_pred_upload["pred_result_val"] = 0
        for i in range(len(Valds_pred_upload)):
            #True Positives
            if Valds_pred_upload["Sample"].iloc[i] == "Pos" and Valds_pred_upload["pred_result"].iloc[i] == "Positive":
                valds_tp += 1
                Valds_pred_upload["pred_result_val"].iloc[i] = "TP"
            #False positives
            if Valds_pred_upload["Sample"].iloc[i] == "Neg" and Valds_pred_upload["pred_result"].iloc[i] == "Positive":
                valds_fp += 1
                Valds_pred_upload["pred_result_val"].iloc[i] = "FP"
            #True Negatives
            if Valds_pred_upload["Sample"].iloc[i] == "Neg" and Valds_pred_upload["pred_result"].iloc[i] == "Negative":
                valds_tn += 1
                Valds_pred_upload["pred_result_val"].iloc[i] = "TN"
            #False negatives
            if Valds_pred_upload["Sample"].iloc[i] == "Pos" and Valds_pred_upload["pred_result"].iloc[i] == "Negative":
                valds_fn += 1
                Valds_pred_upload["pred_result_val"].iloc[i] = "FN"

        valds_fpr = valds_fp/(valds_fp + valds_tn )
        valds_tpr = valds_tp/(valds_tp + valds_fn )

        lista_labels = ["True Positives", "False Negatives", "True Negatives", "False Positives"]
        lista_valores =[valds_tp, valds_fn, valds_tn, valds_fp]

        if pos_lb == 1: 
            df_val_pred_ord = Valds_pred_upload.sort_values(by="Prediction_result", ascending=True)
        if pos_lb == 0:
             df_val_pred_ord = Valds_pred_upload.sort_values(by="Prediction_result", ascending=False)
             
        #st.dataframe(df_val_pred_ord)

        top_10_tp = 0
        for i in range(0, 10):
            if df_val_pred_ord["pred_result_val"].iloc[i] == "TP":
                top_10_tp =+ 1

        top_20_tp = 0
        for i in range(0, 20):
            if df_val_pred_ord["pred_result_val"].iloc[i] == "TP":
                top_20_tp =+ 1

        df_top =df_val_pred_ord.reset_index()
        df_top = df_top.filter(['Sample','Neoantigen','HLA', 'Author', 'Origin', 'Prediction_result', 'pred_result' ], axis=1)
        df_top = df_top.reindex(['Neoantigen','HLA', 'Author', 'Origin','Sample', 'Prediction_result', 'pred_result'], axis=1)
        df_top = df_top.rename(columns={"Sample":"NeoType", "pred_result":"Prediction_Neotype"})
        df_top_10 = df_top.head(10)
        df_top_20 = df_top.head(20)


        
        f1, f2, f3 = st.columns([1,1,3])

        with f1:
            st.metric(label = ("True Positives:"),
                            value = round(valds_tp, 4))
            st.metric(label = ("False Negatives:"),
                            value = round(valds_fn, 4)) 
            st.metric(label = ("True psoitive rate // sensitivity:"),
                            value = round(valds_tpr, 4))  
            st.metric(label = ("Top 10 true positives:"),
                            value = round(top_10_tp, 4))
        with f2:
            st.metric(label = ("True Negatives:"),
                            value = round(valds_tn, 4))
            st.metric(label = ("False Positives:"),
                            value = round(valds_fp, 4))
            st.metric(label = ("False positive rate // (1-specificity):"),
                            value = round(valds_fpr, 4))     
            st.metric(label = ("Top 20 true positives:"),
                            value = round(top_20_tp, 4))           
        # st.write("True Positives: ", valds_tp)
        # st.write("False Positives: ", valds_fp)
        # st.write("True Negatives: ", valds_tn)
        # st.write("False Negatives: ", valds_fn)

        # st.write("True psoitive rate // sensitivity: ", valds_tpr)
        # st.write("False positive rate // (1-specificity) : ", valds_fpr)

    
        # st.write("top 10 true positives: ", top_10_tp)
        # st.write("top 20 true positives: ", top_20_tp)
        with f3:
            with st.expander(label = "Top 10 values", expanded = True):
                st.dataframe(df_top_10, use_container_width=True)
            with st.expander(label = "Top 20 values", expanded = False):
                st.dataframe(df_top_20, use_container_width=True)
            with st.expander(label = "Classification pie chart"):
                pie_class = go.Figure(data=[go.Pie(labels=lista_labels, values=lista_valores)])
                st.plotly_chart(pie_class, use_container_width=True)

        st.header("ICB data set")

        #df_ICB = pd.read_csv("TNB_dataset_MHC_results.csv")

        TNBds_pred_upload["pred_result"] = 0



        if pos_lb == 0:
            for i in range(len(TNBds_pred_upload)):
                if TNBds_pred_upload["Prediction_result"].iloc[i] < dop_thr:
                    TNBds_pred_upload["pred_result"].iloc[i] = "Negative"
                else: 
                    TNBds_pred_upload["pred_result"].iloc[i] = "Positive"
        if pos_lb == 1:
            for i in range(len(TNBds_pred_upload)):
                if TNBds_pred_upload["Prediction_result"].iloc[i] >= dop_thr:
                    TNBds_pred_upload["pred_result"].iloc[i] = "Negative"
                else: 
                    TNBds_pred_upload["pred_result"].iloc[i] = "Positive"
                    

        #st.dataframe(TNBds_pred_upload)


        # Agrupacion y conteo de 
        df_ICB_grouped = TNBds_pred_upload.groupby(['Sample', 'pred_result', 'Cohort', 'Response']).aggregate({'pred_result': 'count'})

        df_ICB_grouped = df_ICB_grouped.add_suffix('_Count').reset_index()

        # st.write("DF group")
        # st.dataframe(df_ICB_grouped)

        lista_cohortes = df_ICB_grouped['Cohort'].unique()
        lista_pvalores = []
        lista_columns= df_ICB_grouped.columns


        df_TNB_group = df_ICB_grouped


        #Lista de pacientes
        lista_pt =  df_TNB_group['Sample'].unique()

        #Df filtrando positivos
        df_TNB_Pos = df_TNB_group[df_TNB_group["pred_result"]=="Positive"]

        #Lista de pacientes positivos 
        lista_pt_pos =  df_TNB_Pos['Sample'].unique()

        lista_sin_pos = []
        for i in lista_pt:
            if not  i in lista_pt_pos:
                lista_sin_pos.append(i)
        # st.write(lista_sin_pos)

        #Me fijo si en el dataframe filtrado por positivos existen todos los pacientes
        #Si hay pacientes que no existen, los busco en la base sin filtrar y agrego una fila

        df_group_added = pd.DataFrame()
        for i in lista_sin_pos:
            for j in range(len(df_TNB_group)):
                if (df_TNB_group["Sample"].iloc[j] == i):
                    new_row = {"Sample": i, "pred_result": "Positive", "Cohort": df_TNB_group["Cohort"].iloc[j], "Response": df_TNB_group["Response"].iloc[j], "pred_result_Count": 0 }
                    new_df = pd.DataFrame([new_row])
                    df_group_added = pd.concat([df_group_added, new_df], ignore_index=True)
                    
        df_TNB_group = pd.concat([df_TNB_group, df_group_added], ignore_index=True)            
        # st.dataframe(df_TNB_group)

        df_TNB_Pos_0 = df_TNB_group[df_TNB_group["pred_result"]=="Positive"]

        df_rn = df_TNB_Pos_0[df_TNB_group["Cohort"]=="Riaz_naive"]
        df_rn_r = df_rn[df_rn["Response"]=="R"]
        # st.write("este")
        # st.dataframe(df_rn_r)


        pos_pred_cohort_responses_dict = {}
        for each in df_TNB_Pos_0['Cohort'].unique():
            pos_pred_cohort_responses_dict[each] = {'R': [], 'NR': [], "p_value": -1}

        for j in range(len(df_TNB_Pos_0)):
            _cohort = df_TNB_Pos_0["Cohort"].iloc[j]
            _response = df_TNB_Pos_0["Response"].iloc[j]
            pos_pred_cohort_responses_dict[_cohort][_response].append(df_TNB_Pos_0["pred_result_Count"].iloc[j])

        #st.write(pos_pred_cohort_responses_dict)



        lista_p_value = []
        lista_cohortes = []

        for each in df_TNB_Pos_0['Cohort'].unique():
            stat, pvalue = stats.mannwhitneyu(pos_pred_cohort_responses_dict[each]["R"], pos_pred_cohort_responses_dict[each]["NR"])
            pos_pred_cohort_responses_dict[each]["p_value"] = pvalue
            lista_cohortes.append(each)
            lista_p_value.append(pvalue)

        lista_05 = []
        for i in range(len(lista_p_value)):
            lista_05.append(0.05)

        # st.write(lista_05)
        # st.write(type(lista_cohortes))
        # st.write(lista_p_value)
        #st.write(pos_pred_cohort_responses_dict) 
        g1, g2 = st.columns(2)
        Scatter_plot = go.Figure()
        Scatter_plot.add_trace(go.Scatter(x=lista_cohortes, y=lista_p_value))
        Scatter_plot.add_shape(type="line",
                                x0=0, y0=0.05, x1=5, y1=0.05,
                                line=dict(color="Black",width=2, dash="dot"))

        with g1: 
            
            st.plotly_chart(Scatter_plot, use_container_width=True)

        with g2:
            for each in df_TNB_Pos_0['Cohort'].unique():
                box_plot = go.Figure()
                box_plot.add_trace(go.Box(y=pos_pred_cohort_responses_dict[each]["R"], 
                                        name = "R",
                                        marker_color = 'indianred'))
                box_plot.add_trace(go.Box(y=pos_pred_cohort_responses_dict[each]["NR"], 
                                        name = "NR",
                                        marker_color = 'lightseagreen'))
                with st.expander(label=each, expanded=False):
                    st.write(each,"'s wilcoxon p-value: ", round(pos_pred_cohort_responses_dict[each]["p_value"], 4))
                    st.plotly_chart(box_plot, use_container_width=True)
