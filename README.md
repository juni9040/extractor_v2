#conda 환경 설치
conda env create -n extractor_v2 --file env.yaml
conda active extractor_v2

python run_extractor.py -t 1 --user LBJ --project 1103 -c 10000

/(user)/(project)/Input/ 폴더 안에 (project)_barcode.csv 파일 넣기
/(user)/(project)/Input/Fastq 폴더 안에 flash 완료한 fastq 파일 넣기

python run_extractor.py -t 1 --user LBJ --project 1103 -c 10000
한 번 더 실행
