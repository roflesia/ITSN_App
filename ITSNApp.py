#pip install streamlit_extras
#pip install streamlit

import streamlit as st
import pandas as pd
import numpy as np
import PIL as pil
from Bio import PDB
import plotly.express as px
# from stmol import showmol
# import py3Dmol
# import requests
# import biotite.structure.io as bsio



#sin copiar 
# def render_mol(pdb):
#     pdbview = py3Dmol.view()
#     pdbview.addModel(pdb,'pdb')
#     pdbview.setStyle({'cartoon':{'color':'spectrum'}})
#     pdbview.setBackgroundColor('white')#('0xeeeeee')
#     pdbview.zoomTo()
#     pdbview.zoom(2, 800)
#     pdbview.spin(True)
#     showmol(pdbview, height = 500,width=800)




#Page set up-----------------------------------------------------------------------------
st.set_page_config(
    page_title="ITSNdb app",
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



#TABS-------------------------------------------------------------------------------------------------------
tab1, tab2, tab3, tab4, tab5 = st.tabs(["HOME", "BD", "GRAPHS", "NeoAg", "pdb"])



#TAB1-------------------------------------------------------------------------------------------------------

with tab1:
    st.title("ITSNdb")
    multi= ''' The ITSNdb is a new neoantigen database with know immunogenic and non immunogenic tumor specific antigenic peptides derived from genomic rearrangements, such as single nucleotide variants (SNVs), that satisfy the following criteria:

    1. The wild type counterpart has been identified in the source protein
    2. The MHC-I presentation has been experimentally validated
    3. The positive or negative immunogenicity has been experimentally validated by, for instance, ELISPOT®
    
    In this sence, all peptides in the database have experimental confirmation of their positive/negative immunogenicity (classified as “Positive” and “Negative” neoantigens respectively) as well as their cell surface presentation.The neoantigens were collected and curated from published articles searched on PubMedTM using “neoantigen'' or “neoepitopes” as keywords. The ITSNdb only includes neoantigens whose inclusion criteria were explicitly described in its reference bibliography.
    '''
    
    #st.markdown(multi)
    
    df1=pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/ITSNdb.csv")
    Columnas=df1.columns
    Types_tumors=df1.value_counts(Columnas[0])
    Authors=df1.value_counts(Columnas[1])
    Neoantigen=df1.value_counts(Columnas[3])
    Wild_type=df1.value_counts(Columnas[4])
    HLA_ser=df1.value_counts(Columnas[5])
    GeneSymbol=df1.value_counts(Columnas[9])

    # st.write("Data content")
    # #Total de datos
    Total_entry= df1.shape[0]
    # st.write("Entry (total):", Total_entry)
    # #cantidad de positivos
    # #cantidad de negativos
    Posneg=df1.value_counts(Columnas[6])  
    # st.write("Positive entries: ", Posneg.iat[0])
    # st.write("Negative entries: ", Posneg.iat[1])
    # #tipos de tumores
    #st.write("Tumor Types: ", Types_tumors.shape[0])
    # #cantidad de archor y non archor
    Ancnon=df1.value_counts(Columnas[8])  
    # st.write("Anchor entries: ", Ancnon.iat[1])
    # st.write("Non Anchor entries: ", Ancnon.iat[0])
    
    #HLA Gen A- Gen B- Gen C

    #Se toman los datos de la columna HLA y se separan en gene alelo y proteina
    HLA_df= df1['HLA']
    HLAGENE=[]
    HLAALLELO=[]
    HLAPROT=[]
    for n in HLA_df:
        #st.write(n)
        Lc=[*n]
        #st.write(Lc)
        #st.write(Lc[4])
        HLAGENE.append(Lc[4])
        HLAALLELO.append(Lc[5]+Lc[6])
        HLAPROT.append(Lc[8]+Lc[9])

    df1["HLA- GENE"]=HLAGENE
    df1["HLA- ALLELE GROUP"]=HLAALLELO
    df1["HLA- PROTEIN"]=HLAPROT
    
    
    Columnas=df1.columns
    Hgene=df1.value_counts(Columnas[11])  
    # st.write("Entries HLA-A: ", Hgene.iat[0])
    # st.write("Entries HLA-B: ", Hgene.iat[1])
    # st.write("Entries HLA-C: ", Hgene.iat[2])

    # style_metric_cards( background_color= "#FFF",
    #                     border_size_px= 1,
    #                     border_color= "#CCC",
    #                     border_radius_px= 5,
    #                     border_left_color = "#9AD8E1",
    #                     box_shadow= True,)


    col1, col2, col3=st.columns([0.7, 0.15, 0.15] )
    # style_metric_cards(background_color= "#FFF",
    #                 border_size_px= 1,
    #                 border_color= "#CCC",
    #                 border_radius_px= 5,
    #                 border_left_color = "#9AD8E1",
    #                 box_shadow= True,)
    with col1:
        st.markdown(multi)
    with col2:
        st.metric(label=("Entry (total):"),
                    value=Total_entry)
        st.metric(label=("Tumor Types: "),
                    value=Types_tumors.shape[0])   
        st.metric(label=("Positive entries: "),
                    value=Posneg.iat[0])
    with col3:
        st.metric(label=("Negative entries: "),
                    value=Posneg.iat[1])
        st.metric(label=("Anchor entries: "),
                    value=Ancnon.iat[1])
        st.metric(label=("Non Anchor entries: "),
                    value=Ancnon.iat[0])
    #Cantidad de Autores diferentes
    
    
    #Graficos
    #Positivos vs Negativos
    #Ubicacion de la mutación
    #Tipos de tumores- cant por tipo de tumor y dif pos y neg
    #cant de archor y no anchor


    #codigo para que los datos se transformen en un link
    # st.dataframe(
    # df1,
    # column_config={
    #     "Paper": st.column_config.LinkColumn()}
    # )
    
    # st.write("grafico con matplot")
    # ploty= Posneg.plot(y='NeoType',
    #                    kind="pie",
    #                    figsize=(2, 2),
    #                    title=("Positive vs Negative"))
    # st.pyplot(fig=ploty.figure,
    #             use_container_width=False)
    
    # st.write("grafico con plotly")
    
    

    #Grafico de torta- positivos vs negativos
    # Posneg=Posneg.reset_index('NeoType')
    # fig = px.pie(Posneg, names='NeoType', values='count', width=500, height=500)
    # st.plotly_chart(fig, theme="streamlit", use_container_width=False)


    #Grafico de torta segun tipos tumorales
    dftumor=df1.groupby(['Tumor']).size().reset_index('Tumor')
    fig3 = px.pie(dftumor, names='Tumor', values=0, width=500, height=500)
    st.plotly_chart(fig3, theme="streamlit", use_container_width=False)

    #Grafico de torta- anchor vs non anchor
    # dfanchor=df1.groupby(['PositionType']).size().reset_index('PositionType')
    # fig4=px.pie(dfanchor, names='PositionType', values=0, width=500, height=500)
    # st.plotly_chart(fig4, theme="streamlit", use_container_width=False)


    Tumor_NeoType=df1.groupby(['Tumor', 'NeoType']).size()
    Tumor_NeoType=Tumor_NeoType.reset_index('NeoType')
    Tumor_NeoType=Tumor_NeoType.reset_index('Tumor')
    Tumor_NeoType=Tumor_NeoType.rename(columns={0:'count'})
    fig2 = px.bar(Tumor_NeoType, x="Tumor", y="count", color="NeoType", title="Tumor_neotype", barmode='relative')
    st.plotly_chart(fig2, theme="streamlit", use_container_width=True)

    HLA_NeoType=df1.groupby(['HLA', 'NeoType']).size()
    HLA_NeoType=HLA_NeoType.reset_index('NeoType')
    HLA_NeoType=HLA_NeoType.reset_index('HLA')
    HLA_NeoType=HLA_NeoType.rename(columns={0:'count'})
    fig3 = px.bar(HLA_NeoType, x="HLA", y="count", color="NeoType", title="HLA_neotype", barmode='relative')
    st.plotly_chart(fig3, theme="streamlit", use_container_width=True)
    

    HLAgene_AlleleGroup=df1.groupby(['HLA- GENE', 'HLA- ALLELE GROUP','HLA- PROTEIN']).size()
    HLAgene_AlleleGroup=HLAgene_AlleleGroup.reset_index('HLA- ALLELE GROUP')
    HLAgene_AlleleGroup=HLAgene_AlleleGroup.reset_index('HLA- GENE')
    HLAgene_AlleleGroup=HLAgene_AlleleGroup.rename(columns={0:'count'})
    fig4= px.bar(HLAgene_AlleleGroup, x="HLA- GENE", y="count", color="HLA- ALLELE GROUP", title="HLAGENE_Neotype", barmode='relative')
    st.plotly_chart(fig4, theme="streamlit", use_container_width=True)


    #este es el quehay que trabajar un poco mas:

    HLAgene_AlleleGroupA=df1[df1['HLA- GENE'].str.contains('A')]
    HLAgene_AlleleGroupA=HLAgene_AlleleGroupA.groupby(['HLA- ALLELE GROUP','HLA- PROTEIN']).size()
    HLAgene_AlleleGroupA=HLAgene_AlleleGroupA.reset_index('HLA- ALLELE GROUP')
    HLAgene_AlleleGroupA=HLAgene_AlleleGroupA.reset_index('HLA- PROTEIN')
    HLAgene_AlleleGroupA=HLAgene_AlleleGroupA.rename(columns={0:'count'})
    fig4= px.bar(HLAgene_AlleleGroupA, x="HLA- ALLELE GROUP", y="count", color="HLA- PROTEIN", title="Hola", barmode='relative')
    st.plotly_chart(fig4, theme="streamlit", use_container_width=True)
    
    
    st.dataframe(data=HLAgene_AlleleGroupA)


#TAB2-------------------------------------------------------------------------------------------------------------
with tab2:
    st.write("Here you can visualize the data from de ITSNdb")
   

    st.write(":arrow_left: Also, you can filter all the data by selecting the variables in the sidebar on the left")



    with st.sidebar:
        
        #-------------------------------------------------------------------------------------
        Tumores= st.multiselect(label="Please choose all the tumors types you want to filter", 
                            options=Types_tumors.index[0:200])
        if not Tumores:
            Tumor_types_sel=Types_tumors.index[0:200]
        else:
            Tumor_types_sel=Tumores

        #-------------------------------------------------------------------------------------
        Autores=st.multiselect(label="Please choose all the authors you want to filter for ", 
                            options=Authors.index[0:200])
        if not Autores:
            Authors_sel=Authors.index[0:200]
        else: 
            Authors_sel=Autores

        #-------------------------------------------------------------------------------------
        NeoAg=st.multiselect(label="Neoantigen", 
                            options=Neoantigen.index[0:200])
        
        if not NeoAg:
            NeoAg_sel=Neoantigen.index[0:200]
        else:
            NeoAg_sel=NeoAg

        #-------------------------------------------------------------------------------------
        Wildtype=st.multiselect(label="WT", 
                            options=Wild_type.index[0:200])
        
        if not Wildtype:
            Wildtype_sel=Wild_type.index[0:200]
        else:
            Wildtype_sel=Wildtype

        #-------------------------------------------------------------------------------------      
        Neotype=st.radio(
            "seleccion neotype positivo- negativo- ambos",
            ["Both","Positive", "Negative"], 
            index=0)
            
        
        if Neotype=="Positive":
            Neotype_sel=["Positive"]
        elif Neotype=="Negative":
            Neotype_sel=["Negative"]
        else:
            Neotype_sel=["Positive","Negative"]

        #-------------------------------------------------------------------------------------
        Pos_type=st.radio(
            "Seleccion de tipo de posición",
            ["Both","Anchor", "Non-anchor"],
            index=0)
        
        if Pos_type=="Anchor":
            Pos_type_sel=["Anchor"]
        elif Pos_type=="Non-anchor":
            Pos_type_sel=["Non-anchor"]
        else:
            Pos_type_sel=["Anchor","Non-anchor"]

        #-------------------------------------------------------------------------------------
        Pos_mut=st. slider(label='rango de posicion de mutacion',
                             min_value=1, 
                             max_value=10,
                             value=(1,10),
                             step=1)

        #-------------------------------------------------------------------------------------
        N_Lenght=st. slider(label='Length',
                             min_value=1, 
                             max_value=10,
                             value=(1,10),
                             step=1)


    #Page-------------------------------------------------------------------------------------
    with tab2:
        df2=df1[df1["Author"].isin(Authors_sel) 
            & df1["Tumor"].isin(Tumor_types_sel) 
            & df1["Neoantigen"].isin(NeoAg_sel)
            & df1["WT"].isin(Wildtype_sel)
            & df1["NeoType"].isin(Neotype_sel)
            & df1["PositionType"].isin(Pos_type_sel)
            & (df1["mutPosition"] >=Pos_mut[0])
            & (df1["mutPosition"]<=Pos_mut[1])
            & (df1["Length"] >=N_Lenght[0])
            & (df1["Length"]<=N_Lenght[1])
            ]
        st.dataframe(data=df2)
        st.write("Elimino algunas columnas")
        #df3= df2["Tumor", "Neoantigen", "WT", "GeneSymbol", "HLA", "mutPosition", "PositionType"]

    with tab4:
        #agregar columna de 0
        df2['Select']=0
        #Agrega la columna de seleccion y desabilita la edicion de las otras
        df4=st.data_editor(
                df2,
                column_order=("Select", "Neoantigen", "WT", "HLA", "GeneSymbol" ),
                column_config={
                        "Select": st.column_config.CheckboxColumn(
                        "Select",
                        help="Select the ones you want to see",
                        default=False),
                        "Author": None,
                        "Paper": None,
                        "Length": None, 
                        "NeoType": None,
                        "mutPosition": None,
                        "HLA- GENE": None,
                        "HLA- ALLELE GROUP": None, 
                        "HLA- PROTEIN": None,


                        },
                disabled=["Tumor"], #Agregar el resto de las columnas para que no se puedan editar
                hide_index=True,
        )
        
        #Filtramos los seleccionados
        st.write("Selected")                                                
        df_select=df4[df4["Select"]==1]
        #st.dataframe(data=df_select)


        #Cantidad de elementos a mostrar
        countrows = len(df_select)
        st.write("Cantidad de elementos seleccionados: ", countrows)

        for i in df_select.index:
            neo_s=df_select.loc[i,"Neoantigen"]
            wt_s=df_select.loc[i,"WT"]
            author_s=df_select.loc[i,"Author"]
            paper_s=df_select.loc[i,"Paper"]
            with st.expander(neo_s):
                st.write("Neoantigen:", neo_s)
                st.write("Wildtype:", wt_s)
                st.write("Author:", author_s)
                st.write("Reference Paper: ", paper_s)
                st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Illustration_HLA-A.png/330px-Illustration_HLA-A.png")


with tab5:
    gene1=("GABPA")


#RUTAS INCOMPLETAS

ruta_af_mut=('ARCHIVOSPDB/AF/MUT')


barra_invertida = '/'
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

        
