import csv
import psycopg2
from datetime import datetime

# === CONFIG ===
INPUT_FILE = r'W:\Instruments & Systems\AVAPS\Software\sonde_counter_tail_nc\counter_all_drops_tail.txt'
DB_CONFIG = {
    'dbname': 'dropsonde_db',
    'user': 'postgres',
    'password': 'ncar',
    'host': 'localhost',
    'port': 5432
}

def process_and_update():
    # Connect to database
    conn = psycopg2.connect(**DB_CONFIG)
    cur = conn.cursor()
    #cur.execute("DELETE FROM dropsonde_data") # This line will clear the table before inserting, comment out if you only want to update


    with open(INPUT_FILE, 'r', newline='') as infile:
        reader = csv.reader(infile)
        for row in reader:
            uid = row[0]
            try:
                droptime_str = uid.replace('_', 'T')
                droptime = datetime.strptime(droptime_str, "%Y%m%dT%H%M%S")
            except ValueError:
                continue  # Skip malformed timestamp
            
            operator = row[3]
            serial = row[5]
            lat = float(row[8])
            lon = float(row[9])
            tail = row[10].strip()

            # Insert or ignore if UID already exists
            cur.execute("""
                INSERT INTO dropsonde_data (uid, operator, serial, lat, lon, tail, droptime)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (uid) DO NOTHING
            """, (uid, operator, serial, lat, lon, tail, droptime))

    conn.commit()
    cur.close()
    conn.close()
    print("Database update complete.")

if __name__ == "__main__":
    process_and_update()
