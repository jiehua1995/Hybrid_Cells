"""
Filename: f:\Genome Biology Paper\raw data\DNA-FISH quantification\quantification.py
Path: f:\Genome Biology Paper\raw data\DNA-FISH quantification
Created Date: Monday, January 13th 2025, 12:20:09 am
Modified Date: Monday, January 13th 2025, 12:20:09 am
Author: Jie Hua
Description: This script helps me to count the probe dots in each cell manually and save the data into a CSV file.
Copyright (c) 2025 Jie Hua<Jie.Hua@lmu.de>
"""

import csv

# Prompt user for Sample and Probe
sample = input("Enter sample name: ")
probe = input("Enter probe name: ")

# Create a list to store numeric inputs
data = []

print("Enter numbers (type 'OVER' to finish):")
while True:
    user_input = input("Enter a number: ")
    if user_input.strip().upper() == "OVER":
        break
    try:
        # Try to convert the input to a number and add to the list
        number = float(user_input)
        data.append(number)
    except ValueError:
        print("Invalid input. Please enter a number or 'OVER' to finish.")

# Define output filename
filename = f"{sample}-{probe}.csv"

# Save data to CSV file
with open(filename, mode="w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(["Sample", "Probe", "Number"])
    # Write each row of data
    for number in data:
        writer.writerow([sample, probe, number])

print(f"Data saved to {filename}!")
