using ADAS
build_ADAS_database()
data = retrieve_ADAS_data("C"; year="latest", type="scd", metastable=false)
show_ADAS_data("C")
data = retrieve_ADAS_data("Si"; year="latest", type="acd", metastable=false)