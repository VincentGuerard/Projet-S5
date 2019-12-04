% Detection des violation du maximum de l'angle
violationAngle      = ySystemeNonLineaire(:,9) > deg2rad(5) | ySystemeNonLineaire(:,8) > deg2rad(5);
violationHauteur1   = ySystemeNonLineaire(:,23) < 0.001 | ySystemeNonLineaire(:,24) < 0.001 | ySystemeNonLineaire(:,25) < 0.001;
violationHauteur2   = ySystemeNonLineaire(:,23) > 0.03 | ySystemeNonLineaire(:,24) > 0.03 | ySystemeNonLineaire(:,25) > 0.03;

if violationAngle
    disp('Violation de la contrainte des angles < 5 degres');
elseif violationHauteur1 | violationHauteur2
    disp('Violation des contraintes de hauteur')
else
    disp('Contraintes rencontrees')
end