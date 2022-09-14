#!/bin/bash

# Convert notebook to python executable:
jupyter nbconvert --to script PMIC-Microscale-analysis_template_tslist.ipynb
# Edit the new python program to remove lines with notebook items:
sed -i 's/get_ipython/#get_ipython/g' PMIC-Microscale-analysis_template_tslist.py

if [ $1 = 'calc_only' ]
then
    sed -i "s/print(output_stats)/print(output_stats)\n'''/g" PMIC-Microscale-analysis_template_tslist.py
    echo "'''" >> PMIC-Microscale-analysis_template_tslist.py
fi
