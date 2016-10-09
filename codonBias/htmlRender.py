"""
generates the HTML table for the data acquired
2016
Synchon Mandal
"""
def generateHTMLTable(webpageName, codonData, aminoAcidData):
    with open('%s.html' %(webpageName),'a') as page:
       htmlFramework = """<!DOCTYPE html>
                          <html>
                            <head>
                                <meta charset="UTF-8">
                                <meta name="description" content="table for genome analysis">
                                <title>Genome Analysis</title>
                                <link rel="stylesheet" type="text/css" href="style.css">
                            </head>
                            <body>"""

       page.write(htmlFramework)

       page.write('<div>')

       page.write('<table>')
       page.write('<caption>Codon Analysis</caption>')
       page.write('<tr><th>Codon</th><th>Occurences</th></tr>')
       for data in codonData:
          page.write('<tr><td>' + data[0] + '</td>')
          page.write('<td>' + str(data[1]) + '</td></tr>')
       page.write('</table>')

       page.write('<br>')

       page.write('<table>')
       page.write('<caption>Amino Acid Analysis</caption>')
       page.write('<tr><th>Amino Acid</th><th>Occurences</th></tr>')
       for data in aminoAcidData:
          page.write('<tr><td>' + data[0] + '</td>')
          page.write('<td>' + str(data[1]) + '</td></tr>')
       page.write('</table>')

       page.write('</div>')

       page.write('<script type="text/javascript" src="http://code.jquery.com/jquery.min.js"></script><script type="text/javscript">$(document).ready(function(){$.getJSON("dnaCodonTable.json");});</script></body></html>')
