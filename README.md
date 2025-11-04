#conda 환경 설치
conda env create -n extractor_v2 --file env.yaml

conda active extractor_v2

python run_extractor.py -t 1 --user LBJ --project 1103 -c 10000

/(user)/(project)/Input/ 폴더 안에 (project)_barcode.csv 파일 넣기
/(user)/(project)/Input/Fastq 폴더 안에 flash 완료한 fastq 파일 넣기

python run_extractor.py -t 1 --user LBJ --project 1103 -c 10000 --matched False
한 번 더 실행

barcode 2개가 서로 pair 되어 있지 않고, 두 barcode 사이의 조합을 확인해야 한다면 
command line에 --matched False 추가
