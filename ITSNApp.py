
#-----------------------------------------------------------------------------------------------------------------
#   Librerias utilizadas

import streamlit as st
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots

import roc_utils as ru
from sklearn import metrics
import scipy.stats as stats #test wilcoxon

import plotly.express as px
import plotly.graph_objects as go

from streamlit_extras.metric_cards import style_metric_cards
# import PIL as pil
# import json

# from Bio import PDB
# import plotly.express as px
# import plotly.graph_objects as go
# from plotly import subplots

# from stmol import showmol
# import py3Dmol

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
#   Set Up de la applicación
st.set_page_config(
    page_title="ITSNdb app",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
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

#   TABS
tab1, tab2, tab3, tab4 = st.tabs(["HOME", "DATA BASE", "NEOANTIGEN INFORMATION", "SOFTWARE VALIDATION"])

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
            st.metric(label = ("Total number of entries:"),
                      value = Total_entries)
            st.metric(label = ("Number of positive entries:"),
                      value = Positive_entries)
            st.metric(label = ("Number of anchor entries:"),
                      value = Anchor_entries)
        with col4:
            st.metric(label = ("Tumor types:"),
                      value = Tumor_types)
            st.metric(label = ("Number of negative entries:"),
                      value = Negative_entries)
            st.metric(label = ("Number of non-anchor entries:"),
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
        st.subheader("Data visualization")
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
        
    
        if (st.session_state["contador"]+1) == 1:            
            st.plotly_chart(fig_HLA, theme="streamlit", use_container_width=True)
            st.caption('Fig 1. Neoantigen frequency by HLA allele, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each HLA allele.')
        if (st.session_state["contador"]+1) == 2:
            st.plotly_chart(fig_tumor, theme="streamlit", use_container_width=True)
            st.caption('Fig 2. Pie chart showing the distribution of tumor types in the database.')
        if (st.session_state["contador"]+1) == 3:
            st.plotly_chart(fig_PositionType_NeoType, use_container_width=True)
            st.caption('Fig 3. Bar chart showing the frequency distribution of anchor and non-anchor position types, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each position type.')
        if (st.session_state["contador"]+1) == 4:
            st.plotly_chart(fig_mutPosition_NeoType, use_container_width=True)
            st.caption('Fig 4. Bar chart showing the frequency of mutations at each position, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each position.')
        if (st.session_state["contador"]+1) == 5:
            st.plotly_chart(fig_Length_NeoType, use_container_width=True)
            st.caption('Fig 5. Bar chart showing the frequency of neoantigen lengths, with further breakdown of immunogenic (positive) and non-immunogenic (negative) status for each length.')
        if (st.session_state["contador"]+1) == 6:
            st.plotly_chart(fig_mat_frec_9, use_container_width=True)
            st.caption('Fig 6. FPM matrices of 9-mer neoantigens: A. Immunogenic, B. Non-immunogenic, C. Difference. Each matrix shows the frequency of each amino acid at each position.')
        if (st.session_state["contador"]+1) == 7:   
            st.plotly_chart(fig_mat_frec_10, use_container_width=True)
            st.caption('Fig 7. FPM matrices of 10-mer neoantigens: A. Immunogenic, B. Non-immunogenic, C. Difference. Each matrix shows the frequency of each amino acid at each position.')

        
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
    label_expander = "**Filters**"
    with st.expander(label = label_expander,
                      expanded = True):
        
    
        C1, C2, C3, C4, C5 = st.columns([2,2,2,1,1])

        with C1:
            Neoag_multiselect = st.multiselect(label = "Neoantigens", 
                                        options = Neoag_list)
            
            HLA_gene_multiselect = st.multiselect(label = "HLA gene",
                                            options = HLA_gene_list)
            
            gene_multiselect = st.multiselect (label = "Gene Symbol",
                                        options = Gene_list)

        with C2:
            Wt_multiselect = st.multiselect(label = "WildType",
                                            options = wt_list)
            HLA_alllele_multiselect = st.multiselect(label = "HLA allele group",
                                                options = HLA_Allele_list)
            Tumor_multiselect = st.multiselect(label = "Tumors", 
                                            options = Tumor_list)
        
        with C3:
            mutpos_multiselect = st.multiselect(label = "Mutation position",
                                        options = Mutpos_list)
            HLA_protein_multiselect = st.multiselect(label = "HLA protein",
                                                options = HLA_protein_list)

        with C4:
            Neotype_list=st.radio(label = "Inmunogenecity",
                            options = ["Both","Positive", "Negative"], 
                            index=0)
            Pos_type=st.radio(label = "Position type",
                            options = ["Both","Anchor", "Non-anchor"],
                            index=0)
        with C5:
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

    st.write("   ")
    
    df_filtered["Select"] = 0
    
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
                                  height = 600
                                    ) 
    

with tab3:
    df_Selected = df_filtered[df_filtered["Select"]==1]


    if df_Selected.empty:
        st.caption("Please, select neoantigens in the database section")
    #st.dataframe(df_Selected)
    for i in df_Selected.index:
            
        #EXTRACCIÓN DE LOS DATOS DEL REGISTRO SELECCIONADO            
        neo_s=df_Selected.loc[i,"Neoantigen"]
        wt_s=df_Selected.loc[i,"WT"]
        author_s=df_Selected.loc[i,"Author"]
        paper_s=df_Selected.loc[i,"Paper"]
        gene_s=df_Selected.loc[i,"GeneSymbol"]
        length_s= df_Selected.loc[i, "Length"]

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

        #SE GENERA UN EXPANDER PARA CADA REGISTRO SELECCIONADO
        with st.expander(neo_s):
            
            listas = [list_Neo, list_espacios, list_Wt]
            df_compared_mut = pd.DataFrame(columns= lista_columnas, index=["NeoAntigen", "", "WildType"] , data=listas)
            
            t1,t2= st.columns(2)
            with t1:
                st.metric(label = ("Neoantigen:"),
                        value = neo_s)
                st.metric(label = ("Wildtype:"),
                        value = wt_s)
            with t2:
                st.dataframe(df_compared_mut)
            
                st.metric(label = ("Length:"),
                        value = length_s)

            st.write("Author:", author_s)
            st.write("Reference Paper: ", paper_s)
            #st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Illustration_HLA-A.png/330px-Illustration_HLA-A.png")
            st.write("Gene:", gene_s)

            st.write(list_Neo)

with tab4:
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





    st.header("Software Validation")

    with st.expander("", expanded=True):

        c51, c52 = st.columns([3,5])

        with c51:
            st.image("proc.png")
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

        st.write(correlation)
        
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

<<<<<<< HEAD
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
=======
ruta_af_mut=('ARCHIVOSPDB/AF/MUT')


barra_invertida = '/'
fin_ruta=('.pdb')
>>>>>>> e78b7fed984aa82b8386670720354ce7e314dc41

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

<<<<<<< HEAD
        with g1: 
            st.write("Wilcoxon p-value:")
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
                    st.write(each,"'s wilcoxon p-value: ", pos_pred_cohort_responses_dict[each]["p_value"])
                    st.plotly_chart(box_plot, use_container_width=True)
=======
        
>>>>>>> e78b7fed984aa82b8386670720354ce7e314dc41
