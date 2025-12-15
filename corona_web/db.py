import mysql.connector
import os

def get_connection():
    return mysql.connector.connect(
        host=os.getenv("MYSQL_HOST", "localhost"), # enter your hostname
        user=os.getenv("MYSQL_USER", ""), # enter your username
        password=os.getenv("MYSQL_PASSWORD", ""), # enter your password
        database=os.getenv("MYSQL_DB", "corona_nrdb"),
        port=int(os.getenv("MYSQL_PORT", "3306")), # enter your port number
    )