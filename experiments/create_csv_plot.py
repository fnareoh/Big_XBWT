""" File to extract the number of runs for each methods in each datasets and
output it to csv format (used later to plot graphs).
it takes as parameter the config file listing the dataset and path of the summary
of runs in that dataset. """

import sys

datasets_path = sys.argv[1]
csv_separator = ";"
lines = [["Datasets", "EBWT+RLO", "EBWT+RLO no $", "XBWT"]]


def find_3_last_nb_run(file):
    """ Finds numbers following the string "Number of runs:" in file return the last 3 in a list"""
    l = file.readline()
    res_nb_runs = []
    while l:
        l_split = l.split(":")
        if len(l_split) > 0 and l_split[0] == "Number of runs":
            assert len(l_split) > 1
            res_nb_runs.append(int(l_split[1]))
        l = file.readline()
    return res_nb_runs[-3:]


config_file = open(datasets_path, "r")
list_dataset = [l.split(csv_separator) for l in config_file.read().splitlines()]
config_file.close()

for name, name_file in list_dataset[1:]:
    file = open(name_file, "r")
    lines.append([name] + find_3_last_nb_run(file))
    file.close()

# OUTPUT to csv
print(*lines, sep="\n")
output_file = open("res_" + datasets_path, "w")
for l in lines:
    for i in range(len(l)):
        output_file.write(str(l[i]))
        if i != len(l) - 1:
            output_file.write(csv_separator)
        else:
            output_file.write("\n")
output_file.close()
