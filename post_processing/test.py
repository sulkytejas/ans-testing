import sys
import requests


response = requests.get("http://localhost:3060/api/drive/complete/")
print(response)
