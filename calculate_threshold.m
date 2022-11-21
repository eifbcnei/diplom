function[threshold] = calculate_threshold(array)
    testarray = real(array);
    pks = findpeaks(testarray);
    maxs = sort(pks,'descend');
    max2 = maxs(2);
    threshold = 0.99 * max2;
end