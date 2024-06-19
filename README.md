# DiffAbundanceToolsValidation
Скрипты и пайплайны для валидации тулов дифф экспрессии

Фигма:
https://www.figma.com/file/JkMEbuoJkOFOA1f49f7iya/Untitled?type=whiteboard&node-id=0-1&t=ufQEzK0eGrF50IC5-0

Проверить, что сборщик файлов нормально работает:
python3 assemble_target_files.py

Просто запустить чтобы оно работало:
snakemake -s Snakefile --cores 60

Построить DAG:
snakemake -s Snakefile --dag -n --cores 60 | dot -Tsvg > dag.svg