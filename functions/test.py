lista = [i for i in range(3)]


def test(lista, valor=len(lista)):
	return valor

valor = test(lista)
print(valor)