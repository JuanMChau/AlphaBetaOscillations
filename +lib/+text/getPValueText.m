function pvalue = getPValueText(number)
%GETPVALUETEXT Summary of this function goes here
    if number < 0.0001
        pvalue = 'p<0.0001';
    else
        % Format the number to 4 decimal places without trailing zeros
        formattedValue = sprintf('%.4f', number);
        formattedValue = regexprep(formattedValue, '0+$', '');
        formattedValue = regexprep(formattedValue, '\.$', '');
        pvalue = ['p=' formattedValue];
    end
end