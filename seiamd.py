#!/usr/bin/env python
#
# Copyright (C) 2023 Laura Martinez-Vazquez (lauramv@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import sys, datetime, csv, json
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from collections import OrderedDict

# Global constants

MIN_DATE = datetime.datetime(2020, 1, 26)
MAX_DATE = datetime.datetime(2022, 9, 26)

DATE_FIRST_DEATH = datetime.datetime(2020, 5, 9)

CORRECTION_DATE = datetime.datetime(2021, 10, 29)
CORRECTION_QUANTITY = 4000

LAST_TICK_DATE = MAX_DATE

# Class that holds the parameters of the model

class Parameters:

    def __init__ (self, date, beta_i, lambda_p, alpha, gamma, mu, sigma):
        self.date = date
        self.beta_i = beta_i
        self.lambda_p = lambda_p
        self.alpha = alpha
        self.gamma = gamma
        self.mu = mu
        self.sigma = sigma

    @property
    def beta_a (self):
        # From [Jammi-2020]:
        return .58 * self.beta_i

    def __str__ (self):
        return f'{self.date}: beta_i = {self.beta_i}, beta_a = {self.beta_a}, lambda_p = {self.lambda_p}, alpha = {self.alpha}, gamma = {self.gamma}, mu = {self.mu}, sigma = {self.sigma}'


# Class that holds a list of parameters and provides the latest parameters for a given date 

class ParameterTable:

    def __init__ (self):
        self.parameter_list = []

    def add_parameters (self, parameters):
        self.parameter_list.append(parameters)

    def get_parameters (self, date):
        result = self.parameter_list[0]
        for parameters in sorted(self.parameter_list, key=lambda x: x.date):
            if parameters.date > date:
                break
            result = parameters
        return result

    def __str__ (self):
        return '\n'.join(str(parameters) for parameters in sorted(self.parameter_list, key=lambda x: x.date))


# Class that holds the official data

class Record:

    def __init__ (self, total_population, date, confirmed, discharged=0, deaths=0):
        self.total_population = total_population
        self.date = date
        self.confirmed = confirmed
        self.discharged = discharged
        self.deaths = deaths
        self.vaccinated1 = 0
        self.vaccinated2 = 0
        self.vaccinated3 = 0
        self.vaccinated4 = 0
        self.vaccinated5 = 0
        self.immunized_6m = 0
        self.gamma = None
        self.sigma = None

    def __str__ (self):
        return f'Record({self.date}, {self.confirmed}, {self.discharged}, {self.deaths}, {self.identified}, {self.vaccinated1}, {self.vaccinated2}, {self.vaccinated3}, {self.vaccinated4}, {self.vaccinated5})'

    def __repr__ (self):
        return str(self)
    
    @property
    def identified (self):
        return self.confirmed - self.discharged - self.deaths

    @property
    def asymptomatic (self):
        return self.total_population * .02 - self.identified * .17

    @property
    def immunized_by_vaccine (self):
        v1 = self.vaccinated1 * 0.52
        v2 = self.vaccinated2 * (0.95 - 0.52)
        v3 = self.vaccinated3 * 0.95
        v4 = self.vaccinated4 * 0.95
        v5 = self.vaccinated5 * 0.95
        return v1 + v2 + v3 + v4 + v5

    @property
    def immunized (self):
        return self.gamma * self.identified + self.gamma * self.asymptomatic + self.immunized_by_vaccine - self.sigma * self.immunized_6m


# Read the data from the official CSV and generate a dictionary of records

def load_records (n, fname_confirmed, fname_discharged, fname_deaths, col_confirmed=14, col_discharged=41, col_deaths=14):
    drecords = OrderedDict()
    total_confirmed = 0
    with open(fname_confirmed) as f:
        for row in csv.reader(f):
            if 'Date' in row[0]:
                continue
            for i in range(len(row)):
                if row[i] == '':
                    row[i] = '0'
            year, month, day = row[0].split('/')
            date = datetime.datetime(int(year), int(month), int(day))
            total_confirmed += int(row[col_confirmed])
            if date == CORRECTION_DATE:
                total_confirmed += CORRECTION_QUANTITY
            if date < MIN_DATE or date > MAX_DATE:
                continue
            drecords[date] = Record(n, date, total_confirmed)
    with open(fname_discharged) as f:
        for row in csv.reader(f):
            if 'Date' in row[0]:
                continue
            for i in range(len(row)):
                if row[i] == '':
                    row[i] = '0'
            year, month, day = row[0].split('/')
            date = datetime.datetime(int(year), int(month), int(day))
            if date not in drecords:
                continue
            drecords[date].discharged = int(row[col_discharged])
    with open(fname_deaths) as f:
        for row in csv.reader(f):
            if 'Date' in row[0]:
                continue
            for i in range(len(row)):
                if row[i] == '':
                    row[i] = '0'
            year, month, day = row[0].split('/')
            date = datetime.datetime(int(year), int(month), int(day))
            if date not in drecords:
                continue
            drecords[date].deaths = int(row[col_deaths])
    return drecords


# Read the official vaccination data from the NDJSON and complete the record objects

def load_vaccinated (drecords, fname, prefecture='13', displacement1=14, displacement2=7):
    with open(fname) as f:
        res = []
        for line in f.readlines():
            res.append(json.loads(line))
    for r in res:
        if r['prefecture'] != prefecture:
            continue
        if r['status'] == 1:
            date = datetime.datetime.strptime(r['date'], '%Y-%m-%d') + datetime.timedelta(days=displacement1)
            if date not in drecords:
                print('load_vaccinated: Date not recorded:', date)
                continue
            drecords[date].vaccinated1 += r['count']
        elif r['status'] > 1:
            date = datetime.datetime.strptime(r['date'], '%Y-%m-%d') + datetime.timedelta(days=displacement2)
            if date not in drecords:
                print('load_vaccinated: Date not recorded:', date)
                continue
            attr = {
                2 : 'vaccinated2',
                3 : 'vaccinated3',
                4 : 'vaccinated4',
                5 : 'vaccinated5'
            }
            drecords[date].__dict__[attr[r['status']]] += r['count']
        else:
            print('Wrong status:', r['status'])


# Load the parameters for each period

def load_parameters (fname):
    ptable = ParameterTable()
    with open(fname) as f:
        for row in csv.reader(f):
            if row[0].startswith('year') or row[0].startswith('#'):
                continue
            for i in range(len(row)):
                if row[i] == '' or row[i] == 'None':
                    row[i] = None
            ptable.add_parameters(Parameters(
                datetime.datetime(int(row[0]), int(row[1]), int(row[2])),
                beta_i=(None if row[3] is None else float(row[3])),
                lambda_p=float(row[4]),
                alpha=float(row[5]),
                gamma=float(row[6]),
                mu=float(row[7]),
                sigma=float(row[8])
            ))
    return ptable


# Assign the value of sigma and gamma

def assign_gamma_sigma (drecords, ptable):
    for record in drecords.values():
        record.gamma = ptable.get_parameters(record.date).gamma
        record.sigma = ptable.get_parameters(record.date).sigma


# Calculate the immunized after six months and complete the records

def fill_immunized_6m (drecords):
    for record in drecords.values():
        tdate = record.date - datetime.timedelta(days=6*30)
        if tdate in drecords:
            record.immunized_6m = drecords[tdate].immunized
        else:
            print('fill_immunized_6m: Date not recorded:', tdate)


# Funtion that receives the values of the comparments of the model and the paramenters, and returns the updated values of the compartments of the model

def model_get_next_values (s_t, e_t, i_t, a_t, m_t, d_t, n_t, m_p6m, p_phi, p):
    s_n = s_t - p.beta_i * s_t / n_t * i_t - p.beta_a * s_t / n_t * a_t - p_phi + p.sigma * m_p6m
    e_n = e_t + p.beta_i * s_t / n_t * i_t + p.beta_a * s_t / n_t * a_t - p.lambda_p * e_t
    i_n = i_t + p.alpha * p.lambda_p * e_t - p.gamma * i_t - p.mu * i_t
    a_n = a_t + (1 - p.alpha) * p.lambda_p * e_t - p.gamma * a_t
    m_n = m_t + p.gamma * i_t + p.gamma * a_t + p_phi - p.sigma * m_p6m
    d_n = d_t + p.mu * i_t
    n_n = s_n + e_n + i_n + a_n + m_n + d_n
    return (s_n, e_n, i_n, a_n, m_n, d_n, n_n)


# Generator that yields values of the compartment i, obtained by applying the model to all dates of the record list provided but using only the input data for the first calculation

def model_get_i_projection (n, records, ptable):
    s_t, e_t, i_t, a_t, m_t, d_t, n_t = (n - records[0].identified, 0, records[0].identified, 0, 0, 0, n)
    yield i_t
    for r in records[1:]:
        s_t, e_t, i_t, a_t, m_t, d_t, n_t = model_get_next_values(
            s_t, e_t, i_t, a_t, m_t, d_t, n_t,
            r.immunized_6m, # m_p6m
            r.immunized_by_vaccine, # p_phi
            ptable.get_parameters(r.date)
        )
        yield i_t


# Function that returns the sum of the absolute error of the estimated values of the compartments i, a, d

def model_get_residuals (parameters, s_t, e_t, i_t, a_t, m_t, d_t, n_t, p, records):
    p.beta_i = parameters[0]
    n_r, i_r, a_r, d_r = 0, 0, 0, 0
    for r in records:
        m_p6m = r.immunized_6m
        p_phi = r.immunized_by_vaccine
        s_t, e_t, i_t, a_t, m_t, d_t, n_t = model_get_next_values(s_t, e_t, i_t, a_t, m_t, d_t, n_t, m_p6m, p_phi, p)
        n_r += abs(n_t - r.total_population)
        i_r += abs(i_t - r.identified)
        a_r += abs(a_t - r.asymptomatic)
        d_r += abs(d_t - r.deaths)
    return n_r + i_r + a_r + d_r


# Function that adjusts the parameter and assigns them to the corresponding "parameters" object 

def adjust_parameters (n, records, ptable):
    s_t, e_t, i_t, a_t, m_t, d_t, n_t = (n - records[0].identified, 0, records[0].identified, 0, 0, 0, n)
    current_records = [records[0]]
    p = ptable.get_parameters(records[0].date)
    for r in records[1:]:
        p_new = ptable.get_parameters(r.date)
        if p_new == p:
            current_records.append(r)
            if r.date != records[-1].date:
                continue
        res = scipy.optimize.least_squares(model_get_residuals, .9, args=(s_t, e_t, i_t, a_t, m_t, d_t, n_t, p, current_records), bounds=(0, 10))
        print(current_records[0].date.strftime('%Y-%m-%d'), '->', current_records[-1].date.strftime('%Y-%m-%d'), ':: Adjusting')
        print(res)
        print()
        p.beta_i = res.x[0]
        for cr in current_records:
            s_t, e_t, i_t, a_t, m_t, d_t, n_t = model_get_next_values(
                s_t, e_t, i_t, a_t, m_t, d_t, n_t,
                cr.immunized_6m,
                cr.immunized_by_vaccine,
                p
            )
        current_records = [r]
        p = p_new


# Plot the recorded and calculated values

def plot_comparison (n, records, ptable):
    x = [r.date for r in records]
    y0 = [r.identified for r in records]
    y1 = list(model_get_i_projection(n, records, ptable))
    plt.figure(figsize=(11,5), dpi=300)
    plt.plot(x, y0, label='Official')
    plt.plot(x, y1, label='Calculated')
    plt.xticks([p.date for p in ptable.parameter_list] + [LAST_TICK_DATE], rotation='vertical')
    plt.yticks([0, 50000, 100000, 150000, 200000, 250000, 300000], ['0', '50k', '100k', '150k', '200k', '250k', '300k'])
    plt.xlim(MIN_DATE - datetime.timedelta(days=10), MAX_DATE + datetime.timedelta(days=10))
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.ylabel('Number of persons')
    plt.tight_layout()
    plt.savefig('results-plot-comparison.pdf')
    #plt.show()


def print_records (records):
    print(80 * '-')
    for r in records:
        print(r)
    print(80 * '-')


def print_latex_table (ptable):
    for i in range(len(ptable.parameter_list)):
    	p = ptable.parameter_list[i]
    	if i < len(ptable.parameter_list) - 1:
    	    dend = ptable.parameter_list[i+1].date - datetime.timedelta(days=1)
    	else:
    	    dend = MAX_DATE
    	print(f'({i+1:.0f}) & {p.date:%Y-%m-%d} $\\to$ {dend:%Y-%m-%d} & ${p.beta_i:.2f}$ & ${p.beta_a:.2f}$ \\\\')


def print_calculated_gamma_mu (records_in):
    records = [r for r in records_in if r.date >= DATE_FIRST_DEATH]
    y = [(records[i].discharged - records[i-1].discharged) / records[i].identified for i in range(1, len(records))]
    print('gamma =', sum(y) / len(y))
    y = [(records[i].deaths - records[i-1].deaths) / records[i].identified for i in range(1, len(records))]
    print('mu =', sum(y) / len(y))


# Main fuction

def main ():
    # Set the population (for Tokyo)
    n = int(sys.argv[1])

    # Load the records
    drecords = load_records(n, *sys.argv[2:5])

    # Fill the data of vacccinated
    load_vaccinated(drecords, sys.argv[5])

    # Load parameters table
    ptable = load_parameters(sys.argv[6])

    # Assign gamma and sigma to individual records
    assign_gamma_sigma(drecords, ptable)

    # Calculate the immunized after six months
    fill_immunized_6m(drecords)

    # Get a list of records
    records = list(drecords.values())

    # Print the record list
    print_records(records)

    # Calculate gamma and mu
    print_calculated_gamma_mu(records)

    # Adjust the missing parameter
    adjust_parameters(n, records, ptable)

    # Plot the real and estimated values
    plot_comparison(n, records, ptable)

    # Print the parameters table
    print(ptable)
    
    # Print the parameters table in LaTeX format
    print_latex_table(ptable)


if __name__ == '__main__':
    main()

