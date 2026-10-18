### Worked examples

```mathematica
In[1]:= StringRotateRight["abcdef", 2]
In[2]:= StringRotateRight["abcdef"]  (* one place by default *)
In[3]:= StringRotateRight["abcdef", -2] === StringRotateLeft["abcdef", 2]  (* the two are mirror images *)
In[4]:= StringRotateRight["abc", 100]  (* 100 mod 3 = 1 *)
In[5]:= StringRotateRight["αβγδ", 1]
In[6]:= StringRotateRight[{"abc"}, 1]  (* not threaded over lists: StringRotateRight::string *)
```

### Notes

- `StringRotateRight[s, n]` is `StringRotateLeft[s, -n]`; both share one kernel
  (`rotate` in `src/strings/stringrotate.c`), reducing `n` modulo the character
  count so any machine integer is accepted.
- A non-string first argument emits `::string`, a non-integer shift `::int`.
- Attributes `{Protected}`.
- See also: `StringRotateLeft`, `RotateRight`.
