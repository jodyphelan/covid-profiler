import pyfastx
import argparse
import csv
from tqdm import tqdm
import re
from datetime import datetime

def check_for_disallowed_countries(string):
    for x in ["pangolin","env","bat","canine","mink","tiger","cat","mouse","lion"]:
        if x in string:
            return True
    return False

def date_qc(date):
    re_obj = re.search("([0-9][0-9][0-9][0-9])-([0-9][0-9])-([0-9][0-9])",date)
    if re_obj:
        collection_date = datetime.strptime(date,"%Y-%m-%d")
        if datetime.now() > collection_date:
            return True
    return False


country2country = {
    'Beijing': 'China',
    'Fujian': 'China',
    'Shandong': 'China',
    'Fuyang': 'China',
    'Wuhan': 'China',
    'Guangzhou': 'China',
    'Guangdong': 'China',
    'Hefei': 'China',
    'Shenzhen': 'China',
    'Jingzhou': 'China',
    'Anhui': 'China',
    'Wuhan-Hu-1': 'China',
    'Kanagawa': 'China',
    'Hangzhou': 'China',
    'Shanghai': 'China',
    'Lishui': 'China',
    'Foshan': 'China',
    'Yunnan': 'China',
    'Zhejiang': 'China',
    'NanChang': 'China',
    'Jiangsu': 'China',
    'Jiujiang': 'China',
    'Nanchang': 'China',
    'Shangrao': 'China',
    'Ganzhou': 'China',
    'Xinyu': 'China',
    'Pingxiang': 'China',
    'Jian': 'China',
    'Tianmen': 'China',
    'Brunei': 'China',
    'Tianmen': 'China',
    'Chongqing': 'China',
    'Sichuan': 'China',
    'Jiangxi': 'China',
    'Henan': 'China',
    'NetherlandsL': 'Netherlands',
    'ITALY': 'Italy',
    'Yichun': 'China',
    'Fuzhou': 'China',
    'Yingtan': 'China',
    'SouthAfrica': 'South Africa',
    'Changzhou': 'China',
    'Harbin': 'China',
    'Andalusia': 'Spain',
    'Shaoxing':'China',
    'Bucuresti':'Romania',
    'Romania ':'Romania',
    ' Russia ':'Russia',
    'Bahrein':'Bahrain',
    'SaudiArabia':'Saudi Arabia',
    'Liaoning':'China',
    'Hunan':'China',
    'Gansu':'China',
    'NorthernIreland': 'Northern Ireland',
    'Qingdao':'China',
    'Weifang':'China',
    # '':'China',
    'Crimea':'Ukraine',
}

country2iso_a3 = {
  "Palestine": "PSE",
  "Malta": "MLT",
  "Botswana": "BWA",
  "Zambia": "ZMB",
  "Cuba": "CUB",
  "Gabon": "GAB",
  "Andorra": "AND",
  "Curacao": "CUW",
  "Aruba": "ABW",
  "Gibraltar": "GIB",
  "Madagascar": "MDG",
  "Dominican Republic": "DOM",
  "Suriname": "SUR",
  "Moldova": "MDA",
  "Montenegro": "MNE",
  "Sierra Leone": "SLE",
  "Faroe Islands": "FRO",
  "Ukraine": "UKR",
  "Reunion": "REU",
  "Guatemala": "GTM",
  "Belize": "BLZ",
  "North Macedonia": "MKD",
  "Mali": "MLI",
  "Bulgaria": "BGR",
  "Bahrain": "BHR",
  "Romania ": "ROU",
  "Mongolia": "MNG",
  "Cyprus": "CYP",
  "Benin": "BEN",
  "Venezuela": "VEN",
  "Bosnia and Herzegovina": "BIH",
  "Timor-Leste": "TLS",
  "Jamaica": "JAM",
  "Morocco": "MAR",
  "Oman": "OMN",
  "Kenya": "KEN",
  "Tunisia": "TUN",
  "Iraq": "IRQ",
  "JAM": "JAM",
  "Myanmar": "MMR",
  "Guam": "GUM",
  "Uzbekistan": "UZB",
  "Lebanon": "LBN",
  "Bangladesh": "BGD",
  "Romania": "ROU",
  "Serbia": "SRB",
  "Uganda": "UGA",
  "Scotland": "GBR",
  "Northern Ireland": "GBR",
  "DRC": "COD",
  "Kazakhstan": "KAZ",
  "Sri Lanka": "LKA",
  "Jordan": "JOR",
  "Qatar": "QAT",
  "Croatia": "HRV",
  "Uruguay": "URY",
  "Egypt": "EGY",
  "Gambia": "GMB",
  "Costa Rica": "CRI",
  "Puerto Rico": "PRI",
  "Brunei": "BRN",
  "United Arab Emirates": "ARE",
  "Wales": "GBR",
  "USA": "USA",
  "United Kingdom": "GBR",
  "Finland": "FIN",
  "Israel": "ISR",
  "Canada": "CAN",
  "Hong Kong": "HKG",
  "England": "GBR",
  "France": "FRA",
  "China": "CHN",
  "Spain": "ESP",
  "Luxembourg": "LUX",
  "Switzerland": "CHE",
  "Iceland": "ISL",
  "Norway": "NOR",
  "Ireland": "IRL",
  "Australia": "AUS",
  "Iran": "IRN",
  "Russia": "RUS",
  "Portugal": "PRT",
  "Ghana": "GHA",
  "Netherlands": "NLD",
  "Japan": "JPN",
  "Germany": "DEU",
  "Belgium": "BEL",
  "Brazil": "BRA",
  "Taiwan": "TWN",
  "India": "IND",
  "Italy": "ITA",
  "Senegal": "SEN",
  "Malaysia": "MYS",
  "Thailand": "THA",
  "Austria": "AUT",
  "Singapore": "SGP",
  "Estonia": "EST",
  "South Korea": "KOR",
  "Democratic Republic of the Congo": "COD",
  "Indonesia": "IDN",
  "Kuwait": "KWT",
  "Saudi Arabia": "SAU",
  "Philippines": "PHL",
  "Argentina": "ARG",
  "Vietnam": "VNM",
  "Congo": "COG",
  "Georgia": "GEO",
  "Lithuania": "LTU",
  "Pakistan": "PAK",
  "South Africa": "ZAF",
  "Turkey": "TUR",
  "Korea": "KOR",
  "Latvia": "LVA",
  "New Zealand": "NZL",
  "Chile": "CHL",
  "Slovenia": "SVN",
  "Slovakia": "SVK",
  "Denmark": "DNK",
  "Belarus": "BLR",
  "Mexico": "MEX",
  "Algeria": "DZA",
  "Peru": "PER",
  "Cambodia": "KHM",
  "Ecuador": "ECU",
  "Colombia": "COL",
  "Greece": "GRC",
  "Panama": "PAN",
  "Nigeria": "NGA",
  "Hungary": "HUN",
  "Czech Republic": "CZE",
  "Sweden": "SWE",
  "Nepal": "NPL",
  "Poland": "POL"
}


country2continent = {
    'Palestine':'Asia',
    'Malta':'Europe',
    'Botswana':'Africa',
    'Zambia':'Africa',
    'Cuba':'North America',
    'Hunan':'Asia',
    'Gabon':'Africa',
    'Andorra':'Europe',
    'Curacao':'South America',
    'Aruba':'South America',
    'Gibraltar':'Europe',
    'Madagascar':'North America',
    'Dominican Republic':'North America',
    'Liaoning':'Asia',
    'Suriname':'South America',
    'Moldova':'Europe',
    'Montenegro':'Europe',
    'Sierra Leone':'Africa',
    'SaudiArabia':'Asia',
    'Faroe Islands':'Europe',
    'Ukraine':'Europe',
    'Reunion':'Africa',
    'Guatemala':'North America',
    'Belize':'North America',
    'Bahrein':'Asia',
    ' Russia ':'Europe',
    'North Macedonia':'Europe',
    'Mali':'Africa',
    'Bulgaria':'Europe',
    'Bahrain':'Asia',
    'Romania ':'Europe',
    'Bucuresti':'Europe',
    'Mongolia':'Asia',
    'Cyprus':'Europe',
    'Shaoxing':'Asia',
    'Benin':'Africa',
    'Venezuela':'South America',
    'Andalusia':'Europe',
    'Bosnia and Herzegovina':'Europe',
    'Timor-Leste':'Asia',
    'Jamaica': 'North America',
    'Harbin': 'Asia',
    'Changzhou': 'Asia',
    'SouthAfrica': 'Africa',
    'Yingtan': 'Asia',
    'Fuzhou': 'Asia',
    'Yichun': 'Asia',
    'Morocco': 'Africa',
    'Oman': 'Asia',
    'ITALY': 'Europe',
    'Kenya': 'Africa',
    'Tunisia': 'Africa',
    'Iraq': 'Asia',
    'JAM': 'North America',
    'Myanmar': 'Asia',
    'Guam': 'Asia',
    'Uzbekistan': 'Asia',
    'Lebanon': 'Asia',
    'Bangladesh': 'Asia',
    'Romania': 'Asia',
    'Serbia': 'Asia',
    'Uganda': 'Africa',
    'Beijing': 'Asia',
    'Scotland': 'Europe',
    'Fujian': 'Asia',
    'Shandong': 'Asia',
    'Fuyang': 'Asia',
    'Northern Ireland': 'Europe',
    'DRC': 'Africa',
    'Kazakhstan': 'Asia',
    'Wuhan': 'Asia',
    'Guangzhou': 'Asia',
    'Guangdong': 'Asia',
    'Hefei': 'Asia',
    'Shenzhen': 'Asia',
    'Jingzhou': 'Asia',
    'Sri Lanka': 'Asia',
    'Jordan': 'Asia',
    'Anhui': 'Asia',
    'Wuhan-Hu-1': 'Asia',
    'Kanagawa': 'Asia',
    'Qatar': 'Asia',
    'Hangzhou': 'Asia',
    'Croatia': 'Europe',
    'Shanghai': 'Asia',
    'Uruguay': 'South America',
    'Lishui': 'Asia',
    'Egypt': 'Africa',
    'Gambia': 'Africa',
    'Foshan': 'Asia',
    'Yunnan': 'Asia',
    'NewZealand': 'Oceania',
    'Zhejiang': 'Asia',
    'NanChang': 'Asia',
    'Costa Rica': 'North America',
    'Jiangsu': 'Asia',
    'Jiujiang': 'Asia',
    'Nanchang': 'Asia',
    'Shangrao': 'Asia',
    'Ganzhou': 'Asia',
    'Xinyu': 'Asia',
    'Pingxiang': 'Asia',
    'Jian': 'Asia',
    'Puerto Rico': 'North America',
    'Tianmen': 'Asia',
    'Brunei': 'Asia',
    'Tianmen': 'Asia',
    'Chongqing': 'Asia',
    'Sichuan': 'Asia',
    'Jiangxi': 'Asia',
    'Henan': 'Asia',
    'United Arab Emirates': 'Asia',
    'Wales': 'Europe',
    'NetherlandsL': 'Europe',
    'USA': 'North America',
    'United Kingdom': 'Europe',
    'Finland': 'Europe',
    'Israel': 'Asia',
    'Canada': 'North America',
    'Hong Kong': 'Asia',
    'England': 'Europe',
    'France': 'Europe',
    'China': 'Asia',
    'Spain': 'Europe',
    'Luxembourg': 'Europe',
    'Switzerland': 'Europe',
    'Iceland': 'Europe',
    'Norway': 'Europe',
    'Ireland': 'Europe',
    'Australia': 'Oceania',
    'Iran': 'Asia',
    'Russia': 'Europe',
    'Portugal': 'Europe',
    'Ghana': 'Africa',
    'Netherlands': 'Europe',
    'Japan': 'Asia',
    'Germany': 'Europe',
    'Belgium': 'Europe',
    'Brazil': 'South America',
    'Taiwan': 'Asia',
    'India': 'Asia',
    'Italy': 'Europe',
    'Senegal': 'Africa',
    'Malaysia': 'Asia',
    'Thailand': 'Asia',
    'Austria': 'Europe',
    'Singapore': 'Asia',
    'Estonia': 'Europe',
    'South Korea': 'Asia',
    'Democratic Republic of the Congo': 'Africa',
    'Indonesia': 'Asia',
    'Kuwait': 'Asia',
    'Saudi Arabia': 'Asia',
    'Philippines': 'Asia',
    'Argentina': 'South America',
    'Vietnam': 'Asia',
    'Congo': 'Africa',
    'Georgia': 'Asia',
    'Lithuania': 'Europe',
    'Pakistan': 'Asia',
    'South Africa': 'Africa',
    'Turkey': 'Europe',
    'Korea': 'Asia',
    'Latvia': 'Europe',
    'New Zealand': 'Oceania',
    'Chile': 'South America',
    'Slovenia': 'Europe',
    'Slovakia': 'Europe',
    'Denmark': 'Europe',
    'Belarus': 'Europe',
    'Mexico': 'North America',
    'Algeria': 'Africa',
    'Peru': 'South America',
    'Cambodia': 'Asia',
    'Ecuador': 'South America',
    'Colombia': 'South America',
    'Greece': 'Europe',
    'Panama': 'North America',
    'Nigeria': 'Africa',
    'Hungary': 'Europe',
    'Czech Republic': 'Europe',
    'Sweden': 'Europe',
    'Nepal': 'Asia',
    'Poland': 'Europe'
}

acgt = set(["A","C","G","T"])

def main(args):
    with open(args.prefix+".fasta","w") as O:
        with open(args.prefix+".meta.csv","w") as M:
            writer = csv.DictWriter(M,fieldnames=["id","iso_a3","country","continent","date","seqlen","missing_fraction"])
            writer.writeheader()
            for entry in tqdm(pyfastx.Fasta(args.fasta, full_name=True)):
                seqname = entry.name
                # print(seqname)
                if check_for_disallowed_countries(seqname): continue
                meta = seqname.split("|")
                if len(meta)!=3:continue
                if len(meta[0].split("/"))==1: continue

                country = meta[0].split("/")[1]
                country = country2country.get(country,country)
                iso_a3 = country2iso_a3[country2country.get(country,country)]
                if country=="": continue
                continent = country2continent[country]

                date = meta[2]
                if date_qc(date)==False: continue

                if entry.end<args.seqlen: continue

                missing_chars = sum([n for d,n in entry.composition.items() if d.upper() not in acgt])
                if missing_chars/entry.end>args.missing: continue

                seqid = meta[1]

                writer.writerow(
                    {
                        "id":seqid,"country": country,"continent":continent,
                        "iso_a3": iso_a3, "date":date, "seqlen": entry.end,
                        "missing_fraction": missing_chars/entry.end,
                    }
                )

                seq = list(entry.seq)
                for pos in [i for i,n in enumerate(seq) if n not in acgt]:
                    seq[pos] = "N"
                seq = "".join(seq)
                O.write(">%s\n%s\n" % (seqid,seq))



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,help='File with samples')
parser.add_argument('--prefix',help='VCF file',required=True)
parser.add_argument('--missing',default=0.1,type=float,help='Max fraction of Ns')
parser.add_argument('--seqlen',default=29000,type=int,help='Minimum sequence length')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)