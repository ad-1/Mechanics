# ISS Pass

import requests
import time

# latitude and longitude
lat = 52.6832172
lon = -6.693738

# api url
url = f'http://api.open-notify.org/iss-pass.json?lat={lat}&lon={lon}'
response = requests.get(url)
data = response.json()

# print next iss passes
for p in data['response']:
    rt = p['risetime']
    print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime(rt)))
