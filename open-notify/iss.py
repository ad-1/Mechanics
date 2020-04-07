# ISS Position Realtime

import requests
import time

url = 'http://api.open-notify.org/iss-now.json'

for i in range(10):
    response = requests.get(url)
    data = response.json()
    print('\n  longitude', data['iss_position']['longitude'])
    print('  latitude', data['iss_position']['latitude'])
    time.sleep(1)  # api rate limit is 1 second
