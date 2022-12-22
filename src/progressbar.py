def bar_x_out_of_y(x, y, text: str='') -> None:
    maxbars = 20
    nbars = int(x / y * maxbars)
    print('\rProcessing: | ' + '█' * nbars + '-' * (maxbars - nbars) + ' |', end=' ')
    print(f'~ {x}/{y} {text}' + ' '*15, end='')
    if x >= y:
        print('')
