#!/usr/bin/python

""""
    @file biharmonic_example.py

    @brief Compute biharmonic2_example using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys
import subprocess

gismo_path = os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np


def biharmonic_example(path_bvp, path_output):
    # [!Run biharmonic2_example]
    proc = subprocess.Popen(["../build/bin/biharmonic2_example", "-o", path_output, "-x", path_bvp])
    proc.wait()
    # [!Run biharmonic2_example]

''' ####### USER INPUT ####### '''
max_id = 18

path_xml_coll = "../python_examples/results/XmlCollection/XmlCollection_input.xml"
path_xml_out ="results/XmlCollection/XmlCollection_output.xml"
''' ##### USER INPUT END ##### '''

file = gs.io.gsXmlCollection("")
file.load(path_xml_coll)

file_col = gs.io.gsXmlCollection(path_xml_out)
for id in range(max_id):
    path_bvp = file.getString(id)
    print(path_bvp)

    path_output = path_bvp.replace("bvp/","results/")
    path_output = path_output[:path_output.find('.xml')] + "-result.xml"

    biharmonic_example(path_bvp, path_output)
    file_col.addFile(path_output)

file_col.save()
print("Finished!")



