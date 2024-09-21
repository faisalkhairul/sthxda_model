import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

CSV_FILE = 'C:/Users/faisa/.cursor-tutor/crypto_prices.csv'

def load_data(symbol):
    try:
        df = pd.read_csv(CSV_FILE)
        df = df[df['symbol'] == symbol]
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        df.set_index('timestamp', inplace=True)
        df = df.sort_index()
        return df['price']
    except Exception as e:
        print(f"Terjadi kesalahan saat membaca file CSV: {e}")
        return pd.Series()

def smooth_data(prices, window_length=11, polyorder=2):
    return pd.Series(savgol_filter(prices, window_length, polyorder), index=prices.index)

def predict_next_price(prices, window=30):
    if len(prices) < window:
        return None
    
    smoothed_prices = smooth_data(prices)
    last_ma = smoothed_prices.iloc[-window:].mean()
    forecast = [last_ma] * 15
    
    return forecast

def visualize_data_and_prediction(prices, forecast, symbol):
    plt.figure(figsize=(12, 6))
    plt.plot(prices.index, prices.values, label='Harga Aktual', color='blue')
    
    smoothed_prices = smooth_data(prices)
    plt.plot(smoothed_prices.index, smoothed_prices.values, label='Harga Diperhalus', color='green')
    
    last_timestamp = prices.index[-1]
    forecast_timestamps = [last_timestamp + timedelta(minutes=i+1) for i in range(15)]
    plt.plot(forecast_timestamps, forecast, label='Prediksi', color='red', linestyle='--')
    
    plt.title(f'Harga dan Prediksi untuk {symbol}')
    plt.xlabel('Waktu')
    plt.ylabel('Harga (USD)')
    plt.legend()
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Contoh penggunaan
crypto_symbols = ['BNB', 'ETH', 'MATIC', 'FTM']
for crypto_symbol in crypto_symbols:
    print(f"\nAnalisis harga untuk {crypto_symbol}:")
    prices = load_data(crypto_symbol)

    if not prices.empty:
        try:
            forecast = predict_next_price(prices)
            if forecast:
                print("\nPrediksi harga untuk 15 menit ke depan:")
                for i, price in enumerate(forecast):
                    print(f"Minute {i+1}: ${price:.4f}")
                
                visualize_data_and_prediction(prices, forecast, crypto_symbol)
            else:
                print("Tidak cukup data untuk melakukan prediksi")
        except Exception as e:
            print(f"Terjadi kesalahan saat memprediksi harga: {e}")
    else:
        print("Data tidak tersedia")

    # Cetak beberapa statistik data
    print(f"\nStatistik data untuk {crypto_symbol}:")
    print(f"Jumlah data: {len(prices)}")
    if not prices.empty:
        print(f"Range waktu: {prices.index.min()} sampai {prices.index.max()}")
        print(f"Harga terakhir: ${prices.iloc[-1]:.4f}")