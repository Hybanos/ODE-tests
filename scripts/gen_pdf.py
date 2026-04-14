import weasyprint
import os

def scaling():
    out = "<h1>Benchmarks</h1>"

    out += "<img src='./scripts/imgs/raw_time.svg'/>"
    out += "<img src='./scripts/imgs/time_per_step.svg'/>"

    return out

def errors():
    out = "<h1>Errors</h1>"
    imgs = os.listdir("scripts/imgs/")
    imgs.sort()

    for img in imgs:
        if not img.startswith("err_"): continue

        _, n, name = img.split(".")[0].split("_")

        # out += f"<h2>{name}</h2>"
        out += f"<img src='./scripts/imgs/{img}'/>"
    
    return out

out = f"""<html>
<head>
<link rel="stylesheet" href="scripts/haha.css"/>
</head>
<body>
{scaling()}
{errors()}
</body></html>
"""

with open("out.html", "w") as f:
    f.write(out)

css = weasyprint.CSS(string="""
@page {size: A4; margin: 1cm;} 
th, td {border: 1px solid black;}
""")
weasyprint.HTML("./out.html").write_pdf("out.pdf", stylesheets=[css])