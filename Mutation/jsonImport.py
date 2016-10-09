"""
json importer
2016
Synchon Mandal
"""

import json

def openJson(jsonFile):

   with open('%s.json' %(str(jsonFile)),'r') as table:
      data = json.load(table)
      return data
