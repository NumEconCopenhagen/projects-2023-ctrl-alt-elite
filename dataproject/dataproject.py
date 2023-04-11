#Collecting data through Energinet's API
import requests
import time
import pandas as pd

response = requests.get(
    url='https://api.energidataservice.dk/dataset/ProductionConsumptionSettlement?offset=0&start=2022-01-01T00:00&end=2023-01-01T00:00&sort=HourUTC%20DESC&timezone=dk'
)

result = response.json()

for k, v in result.items():
    print(k, v)

    time.sleep(20)

records = result.get('records', [])
                                           
print('records:')
for record in records:
    print(' ', record)

result.keys()
type(result['records'])

# Convert to DataFrame
df = pd.DataFrame(result['records'])

