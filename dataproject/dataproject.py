#Collecting data through Energinet's API
import requests
import time
import pandas as pd

# Send a GET request to the API and importing the data
url='https://api.energidataservice.dk/dataset/ProductionConsumptionSettlement?offset=0&start=2022-01-01T00:00&end=2023-01-01T00:00&sort=HourUTC%20DESC&timezone=dk'
params = {
    'offset': 0,
    'start': '2022-01-01T00:00',
    'end': '2023-01-01T00:00',
    'sort': 'HourUTC DESC',
    'timezone': 'dk'
}
response = requests.get(url, params=params)

# Parse the JSON response
result = response.json()

# Limit the output to a certain number of records
limit = 10
records = result.get('records', [])[:limit]

# Extract keys from the result
result.keys()

# Check the type of records
record_type = type(result['records'])

# Convert to DataFrame
df = pd.DataFrame(result['records'])

