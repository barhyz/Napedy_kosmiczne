function interpolant = GetRpatInterpolant(full_filename)

Temp_Struct = importdata(full_filename);
row_headers = str2double(string(Temp_Struct.rowheaders(2:end))); %A hack because matlab cannot convert cell to array properly
column_headers = Temp_Struct.data(1, :);
data = Temp_Struct.data(2:end, :);

interpolant = griddedInterpolant({row_headers, column_headers}, data, 'linear', 'nearest');

end