### Fixed

- The Species and Polities pages answered 404 on MySQL 8 even after a
  population pass: the population status probe aliased a column as
  `generated`, a reserved word in MySQL 8 (not in MariaDB), so the query
  failed and every population page stayed hidden. The alias is now quoted.
