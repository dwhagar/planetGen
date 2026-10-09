### Fixed
- Python and the browser now round half-way numbers the same way (away from zero, on the shortest decimal), so 9.995 reads "10" and 1.005 reads "1.01" in both; Python used to give "9.99" and "1".
- A negative number that rounds to zero prints "0", not "-0", in both.
- A checkout with `core.autocrlf=true` no longer changes the bytes of the lock files and the word list, so their hashes agree between Windows and Linux.
