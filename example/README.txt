1. В папке output_exmpl находятся примеры вывода, смотреть файл alignment_summary.html

2. В папке input_exmpl находятся входные данные, на которых пример вывода был получен

3. Чтобы запустить quast на этих данных введите в консоль команду:
./quast.py -R ./example/input/P.striptis/reference.zip  ./example/input/P.striptis/2*.zip
Вывод нужно искать в папке: 
./quast_results/latest/contig_alignment_plot/
