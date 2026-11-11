'''
    Filename: generate_abstract_booklet.py
    Notes: Generates abstract booklet from a given excel sheet, for the IAU ECRWG conference; runs in conda env weasy312
    Author : Gemini + Ayan
    Created: 06-11-26
    Example: run generate_abstract_booklet.py
'''

import pandas as pd
import html
import os

# --- 1. Load Data ---
directory_path = "/Users/acharyya/Work/conferences/IAU_ECRWG_Nov2026/"
submission_file = os.path.join(directory_path, "Synergies between Observational and Theoretical Astronomy_ Results, Challenges and Future Prospects (odpovede).xlsx")
schedule_file = os.path.join(directory_path, "IAU November conference talk schedule.xlsx")
html_file = os.path.join(directory_path, "IAU ECRWG Abstract Booklet.html")
pdf_file = os.path.join(directory_path, "IAU ECRWG Abstract Booklet.pdf")

# ---------------generate dataframe from excel sheet-------------------
df_requested_talks = pd.read_excel(schedule_file, sheet_name='All abstracts')
df_requested_talks = df_requested_talks.rename(columns=str.strip) # to strp leading and trailing whitespaces in column names
df_requested_talks = df_requested_talks[['ID', 'Full name', 'Affiliation', 'Presentation title', 'Abstract']]

df_schedule = pd.read_excel(schedule_file, sheet_name='Tentative talk schedule Ayan')
df_schedule = df_schedule.rename(columns=str.strip) # to strp leading and trailing whitespaces in column names
df_talks = df_requested_talks[df_requested_talks['ID'].isin(df_schedule['Sl no.'])] # some abstracts submitted for talks made it to the schedule
df_nontalks = df_requested_talks[~(df_requested_talks['ID'].isin(df_schedule['Sl no.']))] # the rest of the abstracts submitted for talks became posters
df_nontalks = df_nontalks[df_nontalks['ID']!='T117'] # removing the abtract about trojan stuff from the talks-converted-to-poster list

# --- Load poster abstracts-----
df_all_submissions = pd.read_excel(submission_file)
df_all_submissions = df_all_submissions.rename(columns=str.strip) # to strp leading and trailing whitespaces in column names
df_requested_posters = df_all_submissions[df_all_submissions['Presentation type'] == 'Poster']
df_requested_posters = df_requested_posters[['Full name', 'Affiliation', 'Presentation title', 'Abstract']]
df_nontalks = df_nontalks.drop(labels=['ID'], axis=1) # so that now this also has the exact same columns as df_poster
df_posters = pd.concat([df_nontalks, df_requested_posters], ignore_index=True)

# ---------------generate booklet from dataframe-----------------------
title_col = [c for c in df_talks.columns if 'title' in c.lower()][0]

# Designated Keynote Speakers
keynote_names = df_schedule[df_schedule['Notes/chat'].str.contains('Keynote speaker')]['Full name'].values

def get_last_name(name):
    if not isinstance(name, str):
        return ""
    parts = name.strip().split()
    return parts[-1] if parts else ""

def process_entries(df, check_keynote=True):
    keynotes = []
    regular = []
    for idx, row in df.iterrows():
        name = str(row['Full name']).strip()
        affil = str(row['Affiliation']).strip()
        title = str(row[title_col]).strip()
        abstract = str(row['Abstract']).strip()
        
        is_kn = check_keynote and any(kn.lower() in name.lower() for kn in keynote_names)
        
        item = {
            'name': name,
            'affil': affil,
            'title': title,
            'abstract': abstract,
            'last_name': get_last_name(name),
            'is_keynote': is_kn
        }
        
        if is_kn:
            keynotes.append(item)
        else:
            regular.append(item)

    keynotes = sorted(keynotes, key=lambda x: (x['last_name'].lower(), x['name'].lower()))
    regular = sorted(regular, key=lambda x: (x['last_name'].lower(), x['name'].lower()))
    return keynotes + regular

# --- 2. Separate and Sort Presenters ---
talk_entries = process_entries(df_talks, check_keynote=True)
poster_entries = process_entries(df_posters, check_keynote=False)

# --- SOC Members Info ---
soc_members = [
    ("Jaroslav Merc", "Czechia", "Chair"),
    ("Ayan Acharyya", "Italy", "Co-Chair"),
    ("Emmanuel Bempong-Manful", "United Kingdom", "Co-Chair"),
    ("Amit Seta", "Australia", "Member"),
    ("Eduardo Monfardini Penteado", "Brazil", "Member"),
    ("Fatemeh Zahra Majidi", "Italy", "Member"),
    ("Giovanni Sabatini", "Italy", "Member"),
    ("Giulia De Somma", "Italy", "Member"),
    ("Reetika Joshi", "United States", "Member"),
    ("Susmita Das", "India", "Member"),
    ("Youssef Moulane", "Morocco", "Member"),
]

# --- 3. Build HTML Booklet ---
html_content = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>IAU ECRWG Abstract Booklet</title>
<style>
  /* Base Page Layout & Footer Page Numbers */
  @page {
    size: A4;
    margin: 20mm 18mm 20mm 18mm;
    
    @bottom-center {
      content: "Page " counter(page) " of " counter(pages);
      font-family: 'Helvetica Neue', Arial, sans-serif;
      font-size: 8.5pt;
      color: #64748b;
    }
  }

  /* Hide Page Number on the First (Cover) Page */
  @page :first {
    @bottom-center {
      content: none;
    }
  }

  body {
    font-family: 'Helvetica Neue', Arial, sans-serif;
    color: #1e293b;
    line-height: 1.5;
    font-size: 10pt;
  }

  /* Page Break Utilities */
  .page-break {
    page-break-before: always;
  }

/* Cover Page Styling */
  .cover-page {
    box-sizing: border-box;
    height: 245mm;
    padding: 20px 10px 10px 10px;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    align-items: center;
    text-align: center;
    background: linear-gradient(
                  rgba(10, 15, 30, 0.75),
                  rgba(10, 15, 30, 0.75)
                ),
                url('https://assets.science.nasa.gov/content/dam/science/missions/webb/science/2024/10/STScI-01JBFEVR4ZB8FCYT95FKG3HP1Y.mp4/jcr:content/renditions/Poster%20Image.png');
    background-size: cover;
    background-position: center;
    background-repeat: no-repeat;
    color: #ffffff;
  }
  .cover-top {
    width: 100%;
  }
  .cover-center {
    width: 100%;
    margin-top: auto;
    margin-bottom: auto;
    padding: 20px 0;
  }
  .cover-org {
    font-size: 9.5pt;
    font-weight: 700;
    text-transform: uppercase;
    color: #38bdf8; /* Sky blue accent */
    letter-spacing: 1px;
    margin-bottom: 12px;
  }
  .cover-title {
    font-size: 22pt;
    font-weight: 800;
    color: #ffffff;
    line-height: 1.3;
    margin-bottom: 16px;
  }
  .cover-subtitle {
    font-size: 14pt;
    font-weight: 600;
    color: #e2e8f0;
    margin-bottom: 24px;
  }
  .cover-logo {
    max-width: 200px;
    height: auto;
    margin: 0 auto;
  }

  /* Bottom Section & Corner Logos */
  .cover-bottom {
    width: 100%;
  }
  .bottom-logos {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    align-items: center;
    width: 100%;
    margin-bottom: 14px;
  }
  .side-logo {
    max-height: 55px;
    width: auto;
  }
  .cover-details {
    font-size: 10.5pt;
    color: #f8fafc;
    border-top: 1px solid rgba(255, 255, 255, 0.3);
    padding-top: 12px;
    width: 100%;
    text-align: center;
  }

  /* SOC Page Styling */
  .soc-page {
    padding-top: 10px;
  }
  .section-title-main {
    font-size: 18pt;
    font-weight: 800;
    color: #0f172a;
    border-bottom: 2px solid #0284c7;
    padding-bottom: 6px;
    margin-bottom: 20px;
  }
  .soc-list {
    margin-top: 20px;
  }
  .soc-item {
    padding: 8px 0;
    border-bottom: 1px solid #f1f5f9;
    font-size: 10.5pt;
  }
  .soc-name {
    font-weight: 700;
    color: #0f172a;
  }
  .soc-role {
    font-size: 9pt;
    color: #0284c7;
    font-weight: 600;
    margin-left: 6px;
  }
  .soc-affil {
    color: #475569;
    font-style: italic;
  }

  /* Section Header Styling */
  .section-header {
    font-size: 16pt;
    font-weight: 800;
    color: #0f172a;
    background-color: #f1f5f9;
    padding: 8px 12px;
    border-left: 4px solid #0284c7;
    margin-top: 10px;
    margin-bottom: 24px;
    page-break-after: avoid;
  }

  /* Abstract Entry Styling */
  .entry {
    margin-bottom: 22px;
    page-break-inside: avoid;
  }
  .speaker-name {
    font-size: 11pt;
    font-weight: 700;
    color: #0f172a;
  }
  .keynote-tag {
    font-weight: 600;
    font-size: 9.5pt;
    color: #2563eb;
    margin-left: 4px;
  }
  .speaker-inst {
    font-size: 8.5pt;
    font-style: italic;
    color: #475569;
    margin-top: 1px;
    margin-bottom: 6px;
  }
  .talk-title {
    font-size: 10.5pt;
    font-weight: 700;
    color: #1e40af;
    margin: 4px 0 6px 0;
  }
  .abstract-text {
    font-size: 9.5pt;
    color: #334155;
    text-align: justify;
    white-space: pre-line;
    margin: 0;
  }
</style>
</head>
<body>

<!-- PAGE 1: COVER PAGE -->
<div class="cover-page">
  <!-- Top Organization Header -->
  <div class="cover-top">
    <div class="cover-org">International Astronomical Union | Executive Committee Working Group Early-Career Members</div>
  </div>

  <!-- Vertically Centered Title Block -->
  <div class="cover-center">
    <h1 class="cover-title">Synergies between Observational and Theoretical Astronomy:<br>Results, Challenges and Future Prospects</h1>
    <div class="cover-subtitle">Abstract Booklet</div>
  </div>

  <!-- Bottom Section: Corner Logos & Dates -->
  <div class="cover-bottom">
    <div class="bottom-logos">
      <!-- Bottom Left Logo -->
      <img class="side-logo" src="https://www.iau.org/images/IAUDefaults/iau_wb.jpg?version=FC78BB6E" alt="IAU Logo" onerror="this.style.display='none'">
      <!-- Bottom Right Logo -->
      <img class="side-logo" src="https://www.iau.org/common/Uploaded%20files/Early-Career%20Members%20logo.png" alt="IAU ECM WG Logo" onerror="this.style.display='none'">
    </div>
    
    <div class="cover-details">
      <strong>17–19 November 2026</strong> | Fully Virtual Conference
    </div>
  </div>
</div>

<!-- PAGE 2: SOC MEMBERS -->
<div class="page-break soc-page">
  <h2 class="section-title-main">Scientific Organizing Committee (SOC)</h2>
  <div class="soc-list">
"""

# --- SOC SECTION ---
for name, affil, role in soc_members:
    role_badge = f'<span class="soc-role">({role})</span>' if role != "Member" else ''
    html_content += f"""
    <div class="soc-item">
      <span class="soc-name">{html.escape(name)}</span> {role_badge} — <span class="soc-affil">{html.escape(affil)}</span>
    </div>
"""

# --- TALK ABSTRACTS SECTION ---
html_content += """
  </div>
</div>

<!-- PAGE 3+: TALK ABSTRACTS -->
<div class="page-break">
  <div class="section-header">Talk Abstracts</div>
"""

for item in talk_entries:
    kn_badge = '<span class="keynote-tag">(Keynote Presenter)</span>' if item['is_keynote'] else ''
    html_content += f"""
<div class="entry">
  <div class="speaker-name">{html.escape(item['name'])} {kn_badge}</div>
  <div class="speaker-inst">{html.escape(item['affil'])}</div>
  <div class="talk-title">{html.escape(item['title'])}</div>
  <div class="abstract-text">{html.escape(item['abstract'])}</div>
</div>
"""

# --- POSTER ABSTRACTS SECTION ---
html_content += """
</div>

<div class="page-break">
  <div class="section-header">Poster Abstracts</div>
"""

for item in poster_entries:
    html_content += f"""
<div class="entry">
  <div class="speaker-name">{html.escape(item['name'])}</div>
  <div class="speaker-inst">{html.escape(item['affil'])}</div>
  <div class="talk-title">{html.escape(item['title'])}</div>
  <div class="abstract-text">{html.escape(item['abstract'])}</div>
</div>
"""

html_content += "</div>\n</body>\n</html>"

with open(html_file, "w", encoding="utf-8") as f:
    f.write(html_content)

print("HTML Booklet generated: abstract_booklet.html")

# --- 4. Optional: Convert directly to PDF via WeasyPrint ---
try:
    import weasyprint
    weasyprint.HTML(html_file).write_pdf(pdf_file)
    print("PDF successfully generated: abstract_booklet.pdf")
except ImportError:
    print("Tip: Install WeasyPrint (`pip install weasyprint`) to auto-export PDF, or simply open 'abstract_booklet.html' in your browser and click Print -> Save as PDF.")