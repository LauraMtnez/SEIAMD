# SEIAMD
The epidemiological compartmental model SEIAMD (Susceptible-Exposed-Identified-Asymptomatic-iMmunized-Deceased) is an extension of the classical model SIR (Susceptible-Infected-Recovered).

# Invocation for Tokyo
python ./seiamd.py 14039312 newly_confirmed_cases_daily.csv requiring_inpatient_care_etc_daily.csv deaths_cumulative_daily.csv prefecture.ndjson data-20220926-parameters.csv
