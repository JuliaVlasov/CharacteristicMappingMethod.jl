function Mcjet = HMapCompose(params, gszo, Moutjet, gszi, Minjet)

Mcjet = HCompose(params, gszo, Moutjet, gszi, Minjet);

Mcjet = Mcjet + Minjet;

end