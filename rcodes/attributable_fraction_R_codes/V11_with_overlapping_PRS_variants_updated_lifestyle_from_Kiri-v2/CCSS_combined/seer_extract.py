# # !pip install tabula-py
# !pip install pikepdf
# ! pip install requests
# ! pip install PyPDF4



import re
import csv
import PyPDF4

pdf_file = open('sitetype.icdo3.20220429.pdf', 'rb')
pdf_reader = PyPDF4.PdfFileReader(pdf_file)

# Create a regular expression pattern to match the site types
# pattern = r"(\d{4}/\d{1,2})\s+(.+)"
# pattern = r"\b(\d{4}/\d{1,2})\s+([A-Za-z\s-]+)"
pattern = r"(\d{4}/\d{1,2})\s+([^.\n]+)"

# Create a CSV file to write the data to
with open('site_types.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Iterate over each page in the PDF file
    for page_num in range(pdf_reader.getNumPages()):
        page = pdf_reader.getPage(page_num)
        text = page.extractText()
        # print(text)

        matches = re.findall(pattern, text)
        for match in matches:
            print(match)
            writer.writerow(match)


# Close the PDF file and the CSV file
pdf_file.close()
csv_file.close()

