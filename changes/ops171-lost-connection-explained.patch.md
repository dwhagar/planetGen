**Headline:** A statement that loses its MySQL connection now leaves a warning with what the server says about it.

### Changed
- **Lost connections are explained (OPS.171).** When a statement run through the retry helper loses its connection, a warning gives the server's uptime (a small figure means it restarted), its idle timeout and its aborted-client count. The sky unit's own explanation is a warning too. generation-output: unchanged
