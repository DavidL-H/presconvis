###############################################################################
'''
 David Lennox-Hvenekilde
 Modified 29/7/22
 A stripped down wrapper for simple queries to the EMBL clustal omega REST API
 https://www.ebi.ac.uk/Tools/common/tools/help/
'''

import requests

baseUrl = u'https://www.ebi.ac.uk/Tools/services/rest/clustalo'

# Requests function
def restRequest(url):
    result = requests.get(url).content
    return result

# Get xml service parameters
def serviceGetParameters():
    requestUrl = baseUrl + u'/parameters'
    xmlDoc = restRequest(requestUrl)
    return xmlDoc

test = serviceGetParameters()
print(test)


# Post job request

'''
POST /run example

curl -X POST --header 'Content-Type: application/x-www-form-urlencoded' --header 'Accept: text/plain' -d 'email=davidlh3%40gmail.com&title=My%20Job%201&guidetreeout=true&sequence=%3EP12996%0AMAHRPRWTLSQVTELFEKPLLDLLFEAQQVHRQHFDPRQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLESARKAKAAGSTRFCMGAAWKNPHERDMPYLEQMVQGVKAMGLEACMTLGTLSESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVRDAGIKVCSGGIVGLGETVKDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPEEDKDLQLFRKLGLNPQQTAVLAGDNEQQQRLEQALMTPDTDEYYNAAAL%0A%3EA7ZY31%0AMAHRPRWTLSQVTELFEKPLLDLLFEAQQVHRQHFDPRQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLESARKAKAAGSTRFCMGAAWKNPHERDMPYLEQMVQGVKAMGLEACMTLGTLSESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVRDAGIKVCSGGIVGLGETVKDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPEEDKDLQLFRKLGLNPQQTAVLAGDNEQQQRLEQALMTPDTDEYYNAAAL%0A%3EB4TQT9%0AMARHPRWTLSQVTELFEKPLLELLFEAQQIHRQHFDPQQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLDSARKAKNAGSTRFCMGAAWKNPHERDMPYLEQIVQGVKAMGLETCMTLGMLNESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVREAGIKVCSGGIVGLGETVTDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPAEDKDLQLFRKLGLNPQQTRVLAGDNEQQQRLEQTLMTPDTDDYYNAAAL%0A%3EB5F072%0AMARHPRWTLSQVTELFEKPLLELLFEAQQIHRQHFDPQQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLDSARKAKNAGSTRFCMGAAWKNPHERDMPYLEKIVQGVKAMGLETCMTLGMLNESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVREAGIKVCSGGIVGLGETVTDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPAEDKDLQLFRKLGLNPQQTRVLAGDNEQQQRLEQTLMTPDTDDYYNAAAL' 'https://www.ebi.ac.uk/Tools/services/rest/clustalo/run
'''

# POST alignment request #######################
# Wow, this mega string post actually works...
headers = {
    'Content-Type': 'application/x-www-form-urlencoded',
    'Accept': 'text/plain'
}
d = 'email=davidlh3%40gmail.com&title=My%20Job%201&guidetreeout=true&sequence=%3EP12996%0AMAHRPRWTLSQVTELFEKPLLDLLFEAQQVHRQHFDPRQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLESARKAKAAGSTRFCMGAAWKNPHERDMPYLEQMVQGVKAMGLEACMTLGTLSESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVRDAGIKVCSGGIVGLGETVKDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPEEDKDLQLFRKLGLNPQQTAVLAGDNEQQQRLEQALMTPDTDEYYNAAAL%0A%3EA7ZY31%0AMAHRPRWTLSQVTELFEKPLLDLLFEAQQVHRQHFDPRQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLESARKAKAAGSTRFCMGAAWKNPHERDMPYLEQMVQGVKAMGLEACMTLGTLSESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVRDAGIKVCSGGIVGLGETVKDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPEEDKDLQLFRKLGLNPQQTAVLAGDNEQQQRLEQALMTPDTDEYYNAAAL%0A%3EB4TQT9%0AMARHPRWTLSQVTELFEKPLLELLFEAQQIHRQHFDPQQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLDSARKAKNAGSTRFCMGAAWKNPHERDMPYLEQIVQGVKAMGLETCMTLGMLNESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVREAGIKVCSGGIVGLGETVTDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPAEDKDLQLFRKLGLNPQQTRVLAGDNEQQQRLEQTLMTPDTDDYYNAAAL%0A%3EB5F072%0AMARHPRWTLSQVTELFEKPLLELLFEAQQIHRQHFDPQQVQVSTLLSIKTGACPEDCKYCPQSSRYKTGLEAERLMEVEQVLDSARKAKNAGSTRFCMGAAWKNPHERDMPYLEKIVQGVKAMGLETCMTLGMLNESQAQRLANAGLDYYNHNLDTSPEFYGNIITTRTYQERLDTLEKVREAGIKVCSGGIVGLGETVTDRAGLLLQLANLPTPPESVPINMLVKVKGTPLADNDDVDAFDFIRTIAVARIMMPTSYVRLSAGREQMNEQTQAMCFMAGANSIFYGCKLLTTPNPAEDKDLQLFRKLGLNPQQTRVLAGDNEQQQRLEQTLMTPDTDDYYNAAAL'

response = requests.post(baseUrl+"/run", headers=headers, data=d)
response.content

# GET job status ################################
jobid = str(response.content)[2:-1]
jobid = "clustalo-R20220729-150319-0510-47676526-p2m"


headers = {
    'Accept': 'text/plain'
}

response = requests.get(baseUrl+"/status/"+jobid)
response.content


# GET results ################################
jobid = str(response.content)[2:-1]

def clustalo_alignment():

# header needs to be modified depending on result type wanted
# See result types https://www.ebi.ac.uk/Tools/common/tools/help/
headers = {
    'Accept': 'text/x-clustalw-alignment'
}
# https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/clustalo-R20220729-150319-0510-47676526-p2m/aln-clustal_num
response = requests.get(baseUrl+"/result/"+jobid+"/aln-clustal_num")

# NEED TO FIX THIS TO PROPERLY WRITE THE CLUSTAL FILE
with open(jobid+".clustal_num", "w") as clustalo_file:
    clustalo_file.write(str(response.content))






