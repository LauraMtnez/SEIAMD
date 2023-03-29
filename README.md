# SEIAMD
The epidemiological compartmental model SEIAMD (Susceptible-Exposed-Identified-Asymptomatic-iMmunized-Deceased) is an extension of the classical model SIR (Susceptible-Infected-Recovered).

# Invocation for Tokyo
python ./seiamd.py 14039312 newly_confirmed_cases_daily.csv requiring_inpatient_care_etc_daily.csv deaths_cumulative_daily.csv prefecture.ndjson data-20220926-parameters.csv

# Data files
- CSV data files from: Ministry of Health, Labor and Welfare, Japan. Open data. Available online: https://www.mhlw.go.jp/stf/covid-19/open-data.html
- The file "prefecture.ndjson" is too big for upload. You can get it from: https://data.vrs.digital.go.jp/vaccination/opendata/latest/prefecture.ndjson.
