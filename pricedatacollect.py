import requests
import time
from datetime import datetime
import csv

def get_crypto_price(crypto_id):
    url = f"https://api.coingecko.com/api/v3/simple/price?ids={crypto_id}&vs_currencies=usd"
    response = requests.get(url)
    data = response.json()
    return data[crypto_id]['usd']

def collect_data():
    crypto_ids = {
        'BNB': 'binancecoin',
        'ETH': 'ethereum',
        'MATIC': 'matic-network',
        'FTM': 'fantom'
    }
    csv_file = 'crypto_prices.csv'
    
    with open(csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['timestamp', 'symbol', 'price'])
        
        print(f"Mulai mengumpulkan data. File akan disimpan di: {csv_file}")
        
        for _ in range(1440):  # 1440 menit dalam sehari
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            for symbol, crypto_id in crypto_ids.items():
                try:
                    price = get_crypto_price(crypto_id)
                    writer.writerow([timestamp, symbol, price])
                    print(f"{timestamp} - {symbol}: ${price}")
                except Exception as e:
                    print(f"Error saat mengambil data untuk {symbol}: {str(e)}")
                time.sleep(10)  # Tunggu 5 detik antara setiap permintaan
            
            file.flush()  # Memastikan data ditulis ke disk
            print(f"Data telah ditulis ke {csv_file}")
            time.sleep(20)  # Tunggu 30 detik sebelum siklus berikutnya

if __name__ == "__main__":
    collect_data()