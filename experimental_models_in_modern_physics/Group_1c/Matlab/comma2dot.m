function comma2dot(file)

text = fileread(strcat(file,'.asc'));
conv_text = strrep(text, ',', '.');

FID = fopen(strcat(file,'.asc'), 'w');
fwrite(FID, conv_text, 'char');
fclose(FID);

end

