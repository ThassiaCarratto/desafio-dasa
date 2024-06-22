import os
from flask import Flask, request, jsonify
from flask_cors import CORS
import subprocess

app = Flask(__name__)
CORS(app)

# Caminhos para pastas e arquivos
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(BASE_DIR, 'data')
RESULTS_FOLDER = os.path.join(BASE_DIR, 'resultados')
INPUT_VCF = os.path.join(UPLOAD_FOLDER, 'NIST.vcf.gz')
OUTPUT_CSV = os.path.join(RESULTS_FOLDER, 'variantes_anotadas.csv')

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return "Nenhum arquivo foi enviado", 400
    
    file = request.files['file']
    if file.filename == '':
        return "Nenhum arquivo selecionado", 400
    
    filepath = INPUT_VCF
    file.save(filepath)
    
    # Executa o Snakemake com o novo arquivo VCF
    try:
        # Comando para executar o Snakemake
        command = ['snakemake', '-s', 'Snakefile', '--cores', '1']
        
        # Definindo o diretório de trabalho para o diretório do script
        cwd = os.path.dirname(os.path.abspath(__file__))
        
        # Executa o comando usando subprocess
        process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Aguarda a conclusão do processo
        stdout, stderr = process.communicate()
        
        # Verifica se houve algum erro
        if process.returncode != 0:
            return f"Erro ao executar Snakemake: {stderr.decode()}", 500
        
        # Lê o arquivo de saída gerado pelo Snakemake
        with open(OUTPUT_CSV, 'r') as file:
            data = file.read()
        
        return data, 200
    
    except Exception as e:
        return str(e), 500

@app.route('/variantes', methods=['GET'])
def get_variantes():
    try:
        with open(OUTPUT_CSV, 'r') as file:
            data = file.read()
        return data, 200
    except FileNotFoundError:
        return "Arquivo de resultados não encontrado.", 404

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')