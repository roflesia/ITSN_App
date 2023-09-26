import streamlit as st
import pandas as pd
import numpy as np
import PIL as pil

import streamlit as st
import pandas as pd
import numpy as np
import PIL as pil


#Page set up-----------------------------------------------------------------------------

tab1, tab2, tab3 = st.tabs(["HOME", "BD", "GRAPHS"])

with tab1:
    st.write("Introduccion?")
    st.write("Here you can visualize the data from de ITSNdb")
    df1=pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/ITSNdb.csv")

with tab2:
    st.write("Here you can visualize the data from de ITSNdb")
    df1=pd.read_csv("https://raw.githubusercontent.com/elmerfer/ITSNdb/main/data/ITSNdb.csv")

    st.write("Also, you can filter all the data by selecting the variables in the sidebar on the left")

    Columnas=df1.columns
    Types_tumors=df1.value_counts(Columnas[0])
    Authors=df1.value_counts(Columnas[1])
    Neoantigen=df1.value_counts(Columnas[3])
    Wild_type=df1.value_counts(Columnas[4])
    HLA=df1.value_counts(Columnas[5])
    GeneSymbol=df1.value_counts(Columnas[9])

    #COLUMNAS
    Col_tumor, Col_authors=st.columns(2)
      
    with Col_tumor:
        #-------------------------------------------------------------------------------------
        Tumores= st.multiselect(label="Please choose all the tumors types you want to filter", 
                            options=Types_tumors.index[0:200])
        if not Tumores:
            Tumor_types_sel=Types_tumors.index[0:200]
        else:
            Tumor_types_sel=Tumores

    with Col_authors:
        #-------------------------------------------------------------------------------------
        Autores=st.multiselect(label="Please choose all the authors you want to filter for ", 
                            options=Authors.index[0:200])
        if not Autores:
            Authors_sel=Authors.index[0:200]
        else: 
            Authors_sel=Autores


    Col_NeoAg, Col_WT=st.columns(2)

    with Col_NeoAg:
        #-------------------------------------------------------------------------------------
        NeoAg=st.multiselect(label="Neoantigen", 
                            options=Neoantigen.index[0:200])
        
        if not NeoAg:
            NeoAg_sel=Neoantigen.index[0:200]
        else:
            NeoAg_sel=NeoAg
    with Col_WT:
        #-------------------------------------------------------------------------------------
        Wildtype=st.multiselect(label="WT", 
                            options=Wild_type.index[0:200])

    
        if not Wildtype:
            Wildtype_sel=Wild_type.index[0:200]
        else:
            Wildtype_sel=Wildtype

    Col_1, Col_2=st.columns(2)

        #-------------------------------------------------------------------------------------      
    with Col_1:
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
    with Col_2:
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
    with tab1:
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
